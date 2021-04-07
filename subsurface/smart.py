import numpy as np
import warnings

from cm4twc import SubSurfaceComponent


class SMART(SubSurfaceComponent):
    """
    The Soil Moisture Accounting and Routing for Transport [SMART] model
    (`Mockler et al., 2016`_) is an enhancement of the SMARG (Soil
    Moisture Accounting and Routing with Groundwater) lumped, conceptual
    rainfall–runoff model developed at National University of Ireland,
    Galway (`Kachroo, 1992`_), and based on the soil layers concept
    (`O'Connell et al., 1970`_; `Nash and Sutcliffe, 1970`_). Separate
    soil layers were introduced to capture the decline with soil depth
    in the ability of plant roots to extract water for evapotranspiration.
    SMARG was originally developed for flow modelling and forecasting
    and was incorporated into the Galway Real-Time River Flow Forecasting
    System [GFFS] (`Goswami et al., 2005`_). The SMART model reorganised
    and extended SMARG to provide a basis for water quality modelling by
    separating explicitly the important flow pathways in a catchment.

    The sub-surface component of SMART comprises the runoff generation
    and land runoff routing processes. This sub-surface component is
    made up of six soil layers of equal depth and five linear reservoirs.
    The six soil layers are vertically connected to allow for percolation
    and evaporation. The five linear reservoirs represent the different
    pathways for land runoff. Note, the river routing of SMART is not
    included in this component.

    .. _`Mockler et al., 2016`: https://doi.org/10.1016/j.cageo.2015.08.015
    .. _`Kachroo, 1972`: https://doi.org/10.1016/0022-1694(92)90150-T
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
    """

    _inputs_info = {
        'surface_area': {
            'units': 'm2',
            'kind': 'static'
        }
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
                   # component states
                   soil_layers, overland_store, drain_store,
                   inter_store, shallow_gw_store, deep_gw_store,
                   **kwargs):
        # initialise soil layers and linear reservoirs
        soil_layers[-1][:] = 0
        overland_store[-1][:] = 0
        drain_store[-1][:] = 0
        inter_store[-1][:] = 0
        shallow_gw_store[-1][:] = 0
        deep_gw_store[-1][:] = 0

    def run(self,
            # from exchanger
            transpiration, evaporation_soil_surface, evaporation_ponded_water,
            throughfall, snowmelt,
            # component inputs
            surface_area,
            # component parameters
            theta_c, theta_h, theta_d, theta_s, theta_z, theta_sk, theta_fk,
            theta_gk,
            # component states
            soil_layers, overland_store, drain_store,
            inter_store, shallow_gw_store, deep_gw_store,
            # component constants
            **kwargs):

        dt = self.timedelta_in_seconds

        # determine excess rain quantity from snowmelt and throughfall fluxes
        excess_rain = (throughfall + snowmelt) * dt

        # determine unmet ET quantity from ET fluxes
        unmet_peva = (transpiration + evaporation_soil_surface +
                      evaporation_ponded_water) * dt

        # determine limiting conditions
        water_limited = np.isnan(excess_rain)
        energy_limited = ~water_limited

        # initialise current soil layers to their level at previous time step
        soil_layers[0][:] = soil_layers[-1]

        # --------------------------------------------------------------
        # during energy-limited conditions
        # >-------------------------------------------------------------
        # calculate total antecedent soil moisture
        soil_water = np.sum(soil_layers[-1], axis=-1)

        # calculate surface runoff using quick runoff parameter H and
        # relative soil moisture content
        theta_h_prime = theta_h * (soil_water / theta_z)
        overland_flow = theta_h_prime * excess_rain  # excess rainfall
        # contribution to quick surface runoff store
        excess_rain -= overland_flow  # remainder that infiltrates

        # calculate percolation through soil layers
        # (from top layer [1] to bottom layer [6])
        layer_capacity = theta_z / 6.0
        for i in range(6):
            layer_prv = soil_layers[-1][..., i]
            layer_crt = soil_layers[0][..., i]

            # determine space available in layer before reaching full capacity
            space_in_layer = layer_capacity - soil_layers[-1][..., i]
            # turn off warnings because np.nan in comparison will raise one
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                # there is enough space in layer to hold entire excess rain
                layer_crt[excess_rain <= layer_prv] = \
                    (layer_prv[excess_rain <= layer_prv] +
                     excess_rain[excess_rain <= layer_prv])
                excess_rain[excess_rain <= layer_prv] = 0.0
                # there is not enough space in layer to hold entire excess rain
                layer_crt[excess_rain > layer_prv] = \
                    layer_capacity[excess_rain > layer_prv]
                excess_rain[excess_rain > layer_prv] -= \
                    space_in_layer[excess_rain > layer_prv]

        # calculate saturation excess from remaining excess rainfall
        # sat. excess contribution (if not 0) to quick soil matrix runoff store
        drain_flow = theta_d * excess_rain
        # sat. excess contribution (if not 0) to slow soil matrix runoff store
        inter_flow = (1.0 - theta_d) * excess_rain
        # calculate leak from soil layers
        # (i.e. piston flow becoming active during rainfall events)
        theta_s_prime = theta_s * (soil_water / theta_z)

        # calculate soil moisture contributions to runoff stores
        shallow_gw_flow = np.where(energy_limited, 0.0, np.nan)
        deep_gw_flow = np.where(energy_limited, 0.0, np.nan)
        for i in range(6):
            layer_crt = soil_layers[0][..., i]

            # leak to interflow
            # soil moisture outflow reducing exponentially downwards
            leak_interflow = layer_crt * (theta_s_prime ** (i + 1))
            # leaking if enough soil moisture
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                enough_moisture = leak_interflow < layer_crt
            inter_flow[enough_moisture] += leak_interflow[enough_moisture]
            layer_crt[enough_moisture] -= leak_interflow[enough_moisture]

            # leak to shallow groundwater flow
            # soil moisture outflow reducing linearly downwards
            leak_shallow_gw_flow = layer_crt * (theta_s_prime / (i + 1))
            # leaking if enough soil moisture
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                enough_moisture = leak_shallow_gw_flow < layer_crt
            shallow_gw_flow[enough_moisture] += \
                leak_shallow_gw_flow[enough_moisture]
            layer_crt[enough_moisture] -= \
                leak_shallow_gw_flow[enough_moisture]

            # leak to deep groundwater flow
            # soil moisture outflow reducing exponentially upwards
            leak_deep_gw_flow = layer_crt * (theta_s_prime ** (6 - i))
            # leaking if enough soil moisture
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                enough_moisture = leak_deep_gw_flow < layer_crt
            deep_gw_flow[enough_moisture] += leak_deep_gw_flow[enough_moisture]
            layer_crt[enough_moisture] -= leak_deep_gw_flow[enough_moisture]
        # -------------------------------------------------------------<

        # --------------------------------------------------------------
        # during water-limited conditions
        # >-------------------------------------------------------------
        # no soil moisture contribution to runoff stores (replace np.nan by 0)
        overland_flow[water_limited] = 0.0
        drain_flow[water_limited] = 0.0
        inter_flow[water_limited] = 0.0
        shallow_gw_flow[water_limited] = 0.0
        deep_gw_flow[water_limited] = 0.0

        # attempt to satisfy PE from soil layers
        # (from top layer [1] to bottom layer [6])
        for i in range(6):
            layer_prv = soil_layers[-1][..., i]
            layer_crt = soil_layers[0][..., i]

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                enough_moisture = unmet_peva <= layer_prv
                not_enough_moisture = unmet_peva > layer_prv
            # there is enough moisture in layer to satisfy unmet ET
            layer_crt[enough_moisture] = \
                (layer_prv[enough_moisture] -
                 unmet_peva[enough_moisture])
            unmet_peva[enough_moisture] = 0.0
            # there is not enough moisture in layer to satisfy unmet ET
            unmet_peva[not_enough_moisture] = \
                (unmet_peva[not_enough_moisture] -
                 layer_crt[not_enough_moisture]) * theta_c
            layer_crt[not_enough_moisture] = 0.0

        # route runoff
        # overland
        overland_runoff = overland_store[-1] / theta_sk
        overland_store[0][:] = (overland_store[-1] + overland_flow -
                                overland_runoff * dt)
        overland_store[0][overland_store < 0] = 0.0
        # drain
        drain_runoff = drain_store[-1] / theta_sk
        drain_store[0][:] = (drain_store[-1] + drain_flow -
                             drain_runoff * dt)
        drain_store[0][drain_store < 0] = 0.0
        # inter
        inter_runoff = inter_store[-1] / theta_fk
        inter_store[0][:] = (inter_store[-1] + inter_flow -
                             inter_runoff * dt)
        inter_store[inter_store < 0] = 0.0
        # shallow groundwater
        shallow_gw_runoff = shallow_gw_store[-1] / theta_gk
        shallow_gw_store[0][:] = (shallow_gw_store[-1] + shallow_gw_flow -
                                  shallow_gw_runoff * dt)
        shallow_gw_store[0][shallow_gw_store < 0] = 0.0
        # deep groundwater
        deep_gw_runoff = deep_gw_store[-1] / theta_gk
        deep_gw_store[0][:] = (deep_gw_store[-1] + deep_gw_flow -
                               deep_gw_runoff * dt)
        deep_gw_store[deep_gw_store < 0] = 0.0

        return (
            # to exchanger
            {
                'surface_runoff': overland_runoff + drain_runoff + inter_runoff,
                'subsurface_runoff': shallow_gw_runoff + deep_gw_runoff,
                'soil_water_stress': np.sum(soil_layers[0], axis=-1) / theta_z
            },
            # component outputs
            {}
        )

    def finalise(self, **kwargs):
        pass
