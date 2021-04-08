import numpy as np

from cm4twc import OpenWaterComponent


class SMART(OpenWaterComponent):
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

    The open water component of SMART consists in routing the streamflow
    through the river network by means of a linear reservoir.

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
    _outputs_info = {
        'outgoing_water_volume_transport_along_river_channel': {
            'units': 'm3 s-1',
            'description': 'streamflow at outlet'
        }
    }

    def initialise(self,
                   # component states
                   river_store,
                   **kwargs):
        # initialise linear reservoir
        river_store[-1][:] = 0

    def run(self,
            # from exchanger
            surface_runoff, subsurface_runoff, evaporation_openwater,
            # component parameters
            theta_rk,
            # component states
            river_store,
            # component constants
            **kwargs):

        dt = self.timedelta_in_seconds

        # provisionally calculate river flow
        river_flow = river_store[-1] / theta_rk

        # provisionally calculate new river store state
        store = (river_store[-1]
                 + ((surface_runoff + subsurface_runoff) - river_flow) * dt)

        # check whether store has gone negative
        river_flow = np.where(store < 0,
                              # allow max outflow at 95% of what was in store
                              0.95 * ((surface_runoff + subsurface_runoff)
                                      + (river_store[-1] / dt)),
                              river_flow)
        river_store[0][:] = (
                river_store[-1]
                + ((surface_runoff + subsurface_runoff) - river_flow) * dt
        )

        return (
            # to exchanger
            {
                'water_level': river_store
            },
            # component outputs
            {
                'outgoing_water_volume_transport_along_river_channel': river_flow
            }
        )

    def finalise(self, **kwargs):
        pass
