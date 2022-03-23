import numpy as np

import unifhy
from unifhy.settings import dtype_float


class SubSurfaceComponent(unifhy.component.SubSurfaceComponent):
    """The Soil Moisture Accounting and Routing for Transport [SMART] model
    (`Mockler et al., 2016 <https://doi.org/10.1016/j.cageo.2015.08.015>`_)
    is a bucket-type rainfall-runoff model.

    SMART is an enhancement of the SMARG (Soil Moisture Accounting and
    Routing with Groundwater) lumped, conceptual rainfall–runoff model
    developed at National University of Ireland, Galway (`Kachroo, 1992`_),
    and based on the soil layers concept (`O'Connell et al., 1970`_;
    `Nash and Sutcliffe, 1970`_). Separate soil layers were introduced
    to capture the decline with soil depth in the ability of plant roots
    to extract water for evapotranspiration. SMARG was originally developed
    for flow modelling and forecasting and was incorporated into the
    Galway Real-Time River Flow Forecasting System [GFFS]
    (`Goswami et al., 2005`_). The SMART model reorganised and extended
    SMARG to provide a basis for water quality modelling by separating
    explicitly the important flow pathways in a catchment.

    The sub-surface component of SMART comprises the runoff generation
    and land runoff routing processes. This sub-surface component is
    made up of six soil layers of equal depth and five linear reservoirs.
    The six soil layers are vertically connected to allow for percolation
    and evaporation. The five linear reservoirs represent the different
    pathways for land runoff. Note, the river routing of SMART is not
    included in this component.

    .. _`Kachroo, 1992`: https://doi.org/10.1016/0022-1694(92)90150-T
    .. _`O'Connell et al., 1970`: https://doi.org/10.1016/0022-1694(70)90221-0
    .. _`Nash and Sutcliffe, 1970`: https://doi.org/10.1016/0022-1694(70)90255-6
    .. _`Goswami et al., 2005`: https://doi.org/10.5194/hess-9-394-2005

    :contributors: Thibault Hallouin [1,2], Eva Mockler [1,3], Michael Bruen [1]
    :affiliations:
        1. Dooge Centre for Water Resources Research, University College Dublin
        2. Department of Meteorology, University of Reading
        3. Ireland's Environmental Protection Agency
    :licence: GPL-3.0
    :copyright: 2020, University College Dublin
    :codebase: https://github.com/unifhy-org/unifhycontrib-smart
    """

    _inwards = {
        'canopy_liquid_throughfall_and_snow_melt_flux',
        'transpiration_flux_from_root_uptake'
    }
    _outwards = {
        'soil_water_stress_for_transpiration',
        'surface_runoff_flux_delivered_to_rivers',
        'net_groundwater_flux_to_rivers'
    }
    _parameters_info = {
        'theta_c': {
            'description': 'evaporation decay coefficient',
            'units': '1'
        },
        'theta_h': {
            'description': 'quick runoff ratio',
            'units': '1'
        },
        'theta_d': {
            'description': 'drain flow ratio',
            'units': '1'
        },
        'theta_s': {
            'description': 'soil outflow coefficient',
            'units': '1'
        },
        'theta_z': {
            'description': 'effective soil depth',
            'units': 'kg m-2'
        },
        'theta_sk': {
            'description': 'surface reservoir residence time',
            'units': 's'
        },
        'theta_fk': {
            'description': 'interflow reservoir residence time',
            'units': 's'
        },
        'theta_gk': {
            'description': 'groundwater reservoir residence time',
            'units': 's'
        }
    }
    _states_info = {
        'soil_layers': {
            'units': 'kg m-2',
            'divisions': 6
        },
        'overland_store': {
            'units': 'kg m-2'
        },
        'drain_store': {
            'units': 'kg m-2'
        },
        'inter_store': {
            'units': 'kg m-2'
        },
        'shallow_gw_store': {
            'units': 'kg m-2'
        },
        'deep_gw_store': {
            'units': 'kg m-2'
        }
    }

    def initialise(self,
                   # component parameters
                   theta_z,
                   # component states
                   soil_layers, overland_store, drain_store,
                   inter_store, shallow_gw_store, deep_gw_store,
                   **kwargs):

        if not self.initialised_states:
            # initialise soil layers to be half full
            soil_layers.set_timestep(-1, theta_z[..., np.newaxis] / 6 / 2)  # kg m-2
        
    def run(self,
            # from exchanger
            canopy_liquid_throughfall_and_snow_melt_flux,
            transpiration_flux_from_root_uptake,
            # component parameters
            theta_c, theta_h, theta_d, theta_s, theta_z, theta_sk,
            theta_fk, theta_gk,
            # component states
            soil_layers, overland_store, drain_store,
            inter_store, shallow_gw_store, deep_gw_store,
            # component constants
            **kwargs):

        dt = self.timedelta_in_seconds

        # determine excess rain quantity from snowmelt and throughfall fluxes
        excess_rain = canopy_liquid_throughfall_and_snow_melt_flux * dt

        # determine unmet ET quantity from ET fluxes
        unmet_peva = transpiration_flux_from_root_uptake * dt

        # determine limiting conditions
        energy_limited = excess_rain > 0
        water_limited = ~energy_limited

        # initialise current soil layers to their level at previous step
        soil_layers = (
            soil_layers.get_timestep(0),
            soil_layers.get_timestep(-1)
        )
        soil_layers[0][:] = soil_layers[-1]

        # --------------------------------------------------------------
        # under energy-limited conditions
        # >-------------------------------------------------------------

        # calculate total antecedent soil moisture
        soil_water = np.sum(soil_layers[-1], axis=-1)

        # calculate surface runoff using quick runoff parameter H and
        # relative soil moisture content
        theta_h_prime = theta_h * (soil_water / theta_z)
        # excess rainfall contribution to quick surface runoff store
        overland_flow = theta_h_prime * excess_rain
        # remainder that infiltrates
        excess_rain -= overland_flow

        # calculate percolation through soil layers
        # (from top layer [1] to bottom layer [6])
        layer_capacity = theta_z / 6.
        for i in range(6):
            layer_level = soil_layers[0][..., i]

            # determine space in layer before reaching full capacity
            layer_space = layer_capacity - layer_level

            enough_space = excess_rain <= layer_space

            # enough space in layer to hold entire excess rain
            layer_level[:] = np.where(
                energy_limited & enough_space,
                layer_level + excess_rain,
                layer_level
            )
            excess_rain[energy_limited & enough_space] = 0.

            # not enough space in layer to hold entire excess rain
            layer_level[:] = np.where(
                energy_limited & ~enough_space,
                layer_capacity,
                layer_level
            )
            excess_rain = np.where(
                energy_limited & ~enough_space,
                excess_rain - layer_space,
                excess_rain
            )

        # calculate saturation excess from remaining excess rainfall
        # sat. excess contrib. (if not 0) to quick soil runoff store
        drain_flow = theta_d * excess_rain
        # sat. excess contrib. (if not 0) to slow soil runoff store
        inter_flow = (1.0 - theta_d) * excess_rain

        # -------------------------------------------------------------<

        # calculate leak from soil layers
        # (i.e. piston flow becoming active during rainfall events)
        theta_s_prime = theta_s * (soil_water / theta_z)

        # calculate soil moisture contributions to runoff stores
        shallow_gw_flow = np.zeros(self.spaceshape, dtype_float())
        deep_gw_flow = np.zeros(self.spaceshape, dtype_float())

        for i in range(6):
            layer_level = soil_layers[0][..., i]

            # leak to interflow
            leak_interflow = np.where(
                energy_limited,
                # soil moisture outflow reducing exponentially downwards
                layer_level * (theta_s_prime ** (i + 1)),
                # no soil moisture contribution to runoff store
                0
            )
            inter_flow += leak_interflow
            layer_level[:] = layer_level - leak_interflow

            # leak to shallow groundwater flow
            leak_shallow_gw_flow = np.where(
                energy_limited,
                # soil moisture outflow reducing linearly downwards
                layer_level * (theta_s_prime / (i + 1)),
                # no soil moisture contribution to runoff store
                0
            )
            shallow_gw_flow += leak_shallow_gw_flow
            layer_level[:] = layer_level - leak_shallow_gw_flow

            # leak to deep groundwater flow
            leak_deep_gw_flow = np.where(
                energy_limited,
                # soil moisture outflow reducing exponentially upwards
                layer_level * (theta_s_prime ** (6 - i)),
                # no soil moisture contribution to runoff store
                0
            )
            deep_gw_flow += leak_deep_gw_flow
            layer_level[:] = layer_level - leak_deep_gw_flow

        # --------------------------------------------------------------
        # under water-limited conditions
        # >-------------------------------------------------------------

        # attempt to satisfy PE from soil layers
        # (from top layer [1] to bottom layer [6])
        for i in range(6):
            layer_level = soil_layers[0][..., i]

            enough_moisture = unmet_peva <= layer_level

            # enough soil moisture in layer
            layer_level[:] = np.where(
                water_limited & enough_moisture,
                layer_level - unmet_peva,
                layer_level
            )
            unmet_peva[water_limited & enough_moisture] = 0.

            # not enough soil moisture in layer
            layer_level[water_limited & ~enough_moisture] = 0.
            unmet_peva = np.where(
                water_limited & ~enough_moisture,
                theta_c * (unmet_peva - layer_level),
                unmet_peva
            )

        # -------------------------------------------------------------<

        # route runoff

        # overland
        overland_store = (
            overland_store.get_timestep(0),
            overland_store.get_timestep(-1)
        )
        overland_runoff = overland_store[-1] / theta_sk
        overland_store[0][:] = (
            overland_store[-1] + overland_flow - overland_runoff * dt
        )
        overland_store[0][:] *= overland_store[0] > 0

        # drain
        drain_store = (
            drain_store.get_timestep(0),
            drain_store.get_timestep(-1)
        )
        drain_runoff = drain_store[-1] / theta_sk
        drain_store[0][:] = (
            drain_store[-1] + drain_flow - drain_runoff * dt
        )
        drain_store[0][:] *= drain_store[0] > 0

        # inter
        inter_store = (
            inter_store.get_timestep(0),
            inter_store.get_timestep(-1)
        )
        inter_runoff = inter_store[-1] / theta_fk
        inter_store[0][:] = (
            inter_store[-1] + inter_flow - inter_runoff * dt
        )
        inter_store[0][:] *= inter_store[0] > 0

        # shallow groundwater
        shallow_gw_store = (
            shallow_gw_store.get_timestep(0),
            shallow_gw_store.get_timestep(-1)
        )
        shallow_gw_runoff = shallow_gw_store[-1] / theta_gk
        shallow_gw_store[0][:] = (
            shallow_gw_store[-1] + shallow_gw_flow
            - shallow_gw_runoff * dt
        )
        shallow_gw_store[0][:] *= shallow_gw_store[0] > 0

        # deep groundwater
        deep_gw_store = (
            deep_gw_store.get_timestep(0),
            deep_gw_store.get_timestep(-1)
        )
        deep_gw_runoff = deep_gw_store[-1] / theta_gk
        deep_gw_store[0][:] = (
            deep_gw_store[-1] + deep_gw_flow - deep_gw_runoff * dt
        )
        deep_gw_store[0][:] *= deep_gw_store[0] > 0

        return (
            # to exchanger
            {
                'surface_runoff_flux_delivered_to_rivers':
                    overland_runoff + drain_runoff + inter_runoff,
                'net_groundwater_flux_to_rivers':
                    shallow_gw_runoff + deep_gw_runoff,
                'soil_water_stress_for_transpiration':
                    np.sum(soil_layers[0], axis=-1) / theta_z
            },
            # component outputs
            {}
        )

    def finalise(self, **kwargs):
        pass
