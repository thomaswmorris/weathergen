About
==========

``weathergen`` generates time-varying weather profiles using a synthesis of in-situ observations and satellite reanalysis estimates of meteorological parameters, and is able to recreate weather patterns allowing for the simulation frameworks like `maria <https://github.com/thomaswmorris/maria>`_ to simulate climatologically-accurate ground-based astronomical observations.

Methodology
-----------

Atmospheric parameters are simulated using three ingredients:

(1) The probability distributions of each parameter for each time of year and time of day. 
(2) An eigenmodal decomposition of how parameters covary.
(3) Characteristic spectra of how each eigenmode changes over time. 

These are all calculated using climatological data from `ERA5 <https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3803>`_, `GML <https://gml.noaa.gov/obop/>`_, `ESO <https://www.eso.org/sci/facilities>`_, `MKWC <http://mkwc.ifa.hawaii.edu>`_, and `Meteostat <https://meteostat.net/en/>`_. For each site, reanalysis weather parameters are compiled and adjusted (e.g. for diurnal and annual trends) with respect to in-situ meteorological trends. 

The generated profiles depend on reanalysis data which are limited to a temporal resolution of an hour, but other means of adjustment are possible depending on the sophistication of the alternate sources of weather data. For example, the level of turbulent fluctuation in total column water vapor (the largest driver of sub-millimeter atmospheric interference) is not stationary from hour to hour.

Usage
=====

Install the package using pip:

.. code-block:: bash
    
    pip install weathergen
    
In Python, import the package and pass a site and an array of times to the ```generate``` function. For example, to simulate hourly weather for Princeton, New Jersey between June 1st and September 1st we would write 

.. code-block:: python

    import numpy as np
    from datetime import datetime
    import weathergen

    t0 = datetime(2022,6,1).timestamp() # June 1, 2022
    t1 = datetime(2022,9,1).timestamp() # August 1, 2022

    gen_times = np.arange(t0, t1, 3600) # three months of hourly data

    weather = weathergen.generate(site='princeton', time=gen_times)

All available sites are outlined in the dataframe ``weathergen.sites``. Note that the supplied year is arbitrary; the underlying model considers only annual and diurnal climatological variations. The supported sites are listed below, and are also stored in ``weathergen.sites``. Specified times should be supplied in Unix time.

The weather parameters are contained in the attributes of the ``weather`` object, e.g. ``weather.temperature`` or ``weather.pressure``. The values are typically two-dimensional with shape ``(n_height, n_time)``, where ``weather.height`` and ``weather.time`` describe the time and height of each dimension, in Unix time and meters above sea level. Single-level parameters are described by a single number for each time and do not have a height dimension. 
    
Note that the supplied year is arbitrary: the underlying model considers only annual and diurnal climatological variations. The supported sites are listed below, and are also stored in ``weathergen.sites``. Specified times should be supplied in Unix time.

The weather parameters are contained in the attributes of the ``weather`` object, e.g. ``weather.air_temp`` or ``weather.pressure``. The values are typically two-dimensional with shape ``(n_times, n_heights)``, where ``weather.time`` and ``weather.height`` describe the time and height of each dimension, in Unix time and meters above sea level. Single-level parameters are described by a single number for each time and do not have a layer dimension. 
