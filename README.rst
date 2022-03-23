A unifhy-compliant version of the rainfall-runoff model SMART
-------------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/unifhycontrib-smart?style=flat-square&color=00b0f0
   :target: https://pypi.python.org/pypi/unifhycontrib-smart
   :alt: PyPI Version
.. image:: https://img.shields.io/badge/dynamic/json?url=https://zenodo.org/api/records/5780112&label=doi&query=doi&style=flat-square&color=00b0f0
   :target: https://zenodo.org/badge/latestdoi/355583432
   :alt: DOI
.. image:: https://img.shields.io/github/workflow/status/unifhy-org/unifhycontrib-smart/Run%20tests?style=flat-square&label=tests
   :target: https://github.com/unifhy-org/unifhycontrib-smart/actions/workflows/run_tests.yml
   :alt: Tests Status

The Soil Moisture Accounting and Routing for Transport [SMART] model
(`Mockler et al., 2016`_) is a bucket-type rainfall-runoff model.

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

The sub-surface component of SMART comprises the runoff generation
and land runoff routing processes. This sub-surface component is
made up of six soil layers of equal depth and five linear reservoirs.
The six soil layers are vertically connected to allow for percolation
and evaporation. The five linear reservoirs represent the different
pathways for land runoff. Note, the river routing of SMART is not
included in this component.

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


How to install
~~~~~~~~~~~~~~

This package is available on PyPI, so you can simply run:

.. code-block:: bash

   pip install unifhycontrib-smart
