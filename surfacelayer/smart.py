import numpy as np
import warnings

from cm4twc import SurfaceLayerComponent
from cm4twc import dtype_float


class SMART(SurfaceLayerComponent):
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

    The surface layer component of SMART consists in meeting the
    potential evapotranspiration demand either with rainfall under
    energy-limited conditions or with soil moisture under water-limited
    conditions. Note, unlike the original SMART model, this component
    calculates the available soil moisture from the soil water stress
    coefficient provided by the sub-surface component - in the original
    SMART model, the available soil moisture is iteratively depreciated
    with soil layer depth. This unavoidable simplification may
    overestimate the soil moisture available compared to the original
    SMART model.

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
        'precipitation': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
        },
        'potential_evapotranspiration': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
        },
        'surface_area': {
            'units': 'm2',
            'kind': 'static'
        }
    }
    _parameters_info = {
        'theta_t': {
            'description': 'rainfall aerial correction factor',
            'units': '1'
        },
        'theta_z': {
            'description': 'effective soil depth',
            'units': 'kg m-2'
        }
    }

    def initialise(self):
        pass

    def run(self,
            # from exchanger
            soil_water_stress, water_level,
            # component inputs
            precipitation, potential_evapotranspiration,
            surface_area,
            # component parameters
            theta_t, theta_z,
            # component states
            # component constants
            **kwargs):

        # apply parameter T to rainfall data (aerial rainfall correction)
        corrected_rain = precipitation * theta_t
        # determine limiting conditions
        rain_minus_peva = corrected_rain - potential_evapotranspiration
        energy_limited = rain_minus_peva >= 0.0
        water_limited = ~energy_limited

        # --------------------------------------------------------------
        # during energy-limited conditions
        # >-------------------------------------------------------------
        effective_rain = np.where(energy_limited,
                                  rain_minus_peva,
                                  np.nan)
        # -------------------------------------------------------------<

        # --------------------------------------------------------------
        # during water-limited conditions
        # >-------------------------------------------------------------
        # ignore cells where there is rain excess
        unmet_peva = np.where(water_limited,
                              -rain_minus_peva,
                              np.nan)
        # provisionally set contribution as total available moisture
        soil_moisture = soil_water_stress * theta_z
        soil_moisture_contribution = np.where(water_limited,
                                              soil_moisture,
                                              np.nan)
        # turn off warnings because np.nan in comparison will raise one
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            # limit contribution to unmet ET where there is moisture excess
            soil_moisture_contribution = (
                np.where(soil_moisture_contribution >= unmet_peva,
                         unmet_peva,
                         soil_moisture_contribution)
            )
        # -------------------------------------------------------------<

        # calculate actual evapotranspiration
        actual_evapotranspiration = (potential_evapotranspiration +
                                     soil_moisture_contribution)

        return (
            # to exchanger
            {
                'throughfall': effective_rain,
                'snowmelt': np.zeros(self.spaceshape, dtype_float()),
                'transpiration': np.zeros(self.spaceshape, dtype_float()),
                'evaporation_soil_surface': soil_moisture_contribution,
                'evaporation_ponded_water': np.zeros(self.spaceshape, dtype_float()),
                'evaporation_openwater': np.zeros(self.spaceshape, dtype_float())
            },
            # component outputs
            {
                'evapotranspiration': actual_evapotranspiration
            }
        )

    def finalise(self, **kwargs):
        pass
