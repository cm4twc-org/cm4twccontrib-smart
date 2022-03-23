import numpy as np

import unifhy


class SurfaceLayerComponent(unifhy.component.SurfaceLayerComponent):
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

    The surface layer component of SMART consists in meeting the
    potential evapotranspiration demand either with rainfall alone under
    energy-limited conditions or with rainfall and soil moisture under
    water-limited conditions – throughfall is only generated under
    energy-limited conditions. Note, unlike the original SMART model,
    this component calculates the available soil moisture from the soil
    water stress coefficient provided by the sub-surface component – in
    the original SMART model, the available soil moisture is iteratively
    depreciated with soil layer depth. This unavoidable simplification
    may overestimate the soil moisture available compared to the original
    SMART model.

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
        'soil_water_stress_for_transpiration'
    }
    _outwards = {
        'canopy_liquid_throughfall_and_snow_melt_flux',
        'transpiration_flux_from_root_uptake'
    }
    _inputs_info = {
        'rainfall_flux': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
        },
        'potential_water_evapotranspiration_flux': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
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
    _outputs_info = {
        'actual_water_evapotranspiration_flux': {
            'units': 'kg m-2 s-1'
        }
    }

    def initialise(self, **kwargs):
        pass

    def run(self,
            # from exchanger
            soil_water_stress_for_transpiration,
            # component inputs
            rainfall_flux, potential_water_evapotranspiration_flux,
            # component parameters
            theta_t, theta_z,
            # component states
            # component constants
            **kwargs):

        dt = self.timedelta_in_seconds

        # apply parameter T to rainfall data (aerial rainfall correction)
        corrected_rain = rainfall_flux * theta_t

        # determine limiting conditions
        rain_minus_peva = (
            corrected_rain - potential_water_evapotranspiration_flux
        )
        energy_limited = rain_minus_peva > 0.0
        water_limited = ~energy_limited

        # --------------------------------------------------------------
        # under energy-limited conditions
        # >-------------------------------------------------------------

        effective_rain = np.where(
            energy_limited, rain_minus_peva, 0.0
        )

        # -------------------------------------------------------------<

        # --------------------------------------------------------------
        # under water-limited conditions
        # >-------------------------------------------------------------

        # ignore cells where there is rain excess
        unmet_peva = np.where(
            water_limited, -rain_minus_peva, 0.0
        )

        # provisionally set soil evaporation as total available moisture
        soil_water = soil_water_stress_for_transpiration * theta_z
        max_soil_evaporation = np.where(
            water_limited, soil_water / dt, 0.0
        )

        # limit contribution to unmet ET where there is moisture excess
        soil_evaporation = np.where(
            max_soil_evaporation >= unmet_peva,
            unmet_peva,
            max_soil_evaporation
        )

        # -------------------------------------------------------------<

        # calculate actual evapotranspiration
        actual_evapotranspiration = np.where(
            energy_limited,
            potential_water_evapotranspiration_flux,
            corrected_rain + soil_evaporation
        )

        return (
            # to exchanger
            {
                'canopy_liquid_throughfall_and_snow_melt_flux':
                    effective_rain,
                'transpiration_flux_from_root_uptake':
                    soil_evaporation
            },
            # component outputs
            {
                'actual_water_evapotranspiration_flux': 
                    actual_evapotranspiration
            }
        )

    def finalise(self, **kwargs):
        pass
