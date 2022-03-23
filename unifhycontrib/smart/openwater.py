import numpy as np

import unifhy


class OpenWaterComponent(unifhy.component.OpenWaterComponent):
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

    The open water component of SMART consists in routing the streamflow
    through the river network by means of a linear reservoir.

    .. _`Mockler et al., 2016`: https://doi.org/10.1016/j.cageo.2015.08.015
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
        'surface_runoff_flux_delivered_to_rivers',
        'net_groundwater_flux_to_rivers'
    }
    _outwards = {}
    _parameters_info = {
        'theta_rk': {
            'description': 'channel reservoir residence time',
            'units': 's'
        }
    }
    _states_info = {
        'river_store': {
            'units': 'kg m-2'
        }
    }
    _constants_info = {
        'rho_water': {
            'description': 'volumetric mass density of liquid water',
            'units': 'kg m-3',
            'default_value': 1e3
        }
    }
    _outputs_info = {
        'outgoing_water_volume_transport_along_river_channel': {
            'units': 'm3 s-1',
            'description': 'streamflow at outlet'
        }
    }
    _requires_cell_area = True

    def initialise(self,
                   # component states
                   river_store,
                   **kwargs):

        if not self.initialised_states:
            river_store.set_timestep(-1, 0.)

    def run(self,
            # from exchanger
            surface_runoff_flux_delivered_to_rivers, 
            net_groundwater_flux_to_rivers, 
            # component parameters
            theta_rk,
            # component states
            river_store,
            # component constants
            rho_water,
            **kwargs):

        dt = self.timedelta_in_seconds

        total_runoff = (
            surface_runoff_flux_delivered_to_rivers 
            + net_groundwater_flux_to_rivers
        )

        # provisionally calculate river flow
        river_store = (
            river_store.get_timestep(0),
            river_store.get_timestep(-1)
        )
        river_outflow = river_store[-1] / theta_rk

        # provisionally calculate new river store state
        store = river_store[-1] + (total_runoff - river_outflow) * dt

        # check whether store has gone negative
        river_outflow = np.where(
            store < 0,
            # allow max outflow at 95% of what was in store
            0.95 * (total_runoff + river_store[-1] / dt),
            river_outflow
        )
        river_store[0][:] = (
            river_store[-1] + (total_runoff - river_outflow) * dt
        )

        # convert [kg m-2 s-1] to [m3 s-1]
        river_outflow = river_outflow / rho_water * self.spacedomain.cell_area

        return (
            # to exchanger
            {},
            # component outputs
            {
                'outgoing_water_volume_transport_along_river_channel': 
                    river_outflow
            }
        )

    def finalise(self, **kwargs):
        pass
