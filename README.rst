weathergen
==========

``weathergen`` generates time-varying weather profiles using a synthesis of in-situ observations and satellite reanalysis estimates of meteorological parameters, and is able to faithfully recreate weather patterns allowing for the simulation frameworks like `maria <https://github.com/thomaswmorris/maria>`_ to simulate climatologically-accurate ground-based astronomical observations.

Methodology
-----------

``weathergen`` includes climatological data from `ERA5 <https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3803>`_, `GML <https://gml.noaa.gov/obop/>`_, `ESO <https://www.eso.org/sci/facilities>`_, `MKWC <http://mkwc.ifa.hawaii.edu>`_, and `Meteostat <https://meteostat.net/en/>`_. For each site, reanalysis weather parameters are compiled and adjusted (e.g. for diurnal and annual trends) with respect to in-situ meteorological trends. 

After Gaussianizing the atmospheric features, the fluctuations are orthogonalized as to assume each eigenmode as being stationary and statistically independent. We may then use their characteristic spectra to generate novel time series, which can be sent through the aforementioned process in reverse to yield simulated weather data. 

The generated profiles depend on reanalysis data which are limited to a temporal resolution of an hour, but other means of adjustment are possible depending on the sophistication of the alternate sources of weather data. For example, the level of turbulent fluctuation in total column water vapor (the largest driver of sub-millimeter atmospheric interference) is not stationary from hour to hour.

Parameters
----------

``weathergen`` simulates the following weather parameters:

.. list-table:: Single-level parameters
   :widths: 25 50 25
   :header-rows: 1

   * - tag
     - parameter
     - units
   * - total_water_vapor
     - Total water vapor in the atmospheric column.
     - mm of liquid water equivalent
   * - total_liquid_water
     - Total liquid water in the atmospheric column.
     - mm of liquid water equivalent
   * - total_snow_water
     - Total snow water in the atmospheric column.
     - mm of liquid water equivalent
   * - total_ice_water
     - Total ice water in the atmospheric column.
     - mm of liquid water equivalent
   * - total_cloud_cover
     - The proportion of total column covered by clouds 
     - [none]
   * - total_precipitation
     - The rate of precipitation
     - mm/hr
     
.. list-table:: Multi-level parameters
   :widths: 25 50 25
   :header-rows: 1

   * - tag
     - parameter
     - units
   * - air_temp
     - Ambient air temperature 
     - K
   * - pressure
     - Air pressure
     - Pa
   * - wind_north
     - Northward component of mean wind vector
     - m/s
   * - wind_east
     - Eastward component mean wind vector
     - m/s
   * - abs_hum
     - Water vapor content of air
     - kg/m3
   * - ozone
     - Ozone content of air 
     - kg/m3
   * - cloud_cover
     - The proportion of the layer covered by clouds 
     - [none]
   * - divergence
     - Divergence
     - 1/s
   * - vorticity
     - Potential vorticity
     - 1/s

Usage
-----

Install the package using pip:

.. code-block:: bash
    
    pip install weathergen
       
In Python, import the package and pass a site and an array of times to the ``generate`` function. For example, to simulate weather for Princeton, New Jersey between June 1st and September 1st at a resolution of a minute we would write 

.. code-block:: python

    import numpy as np
    from datetime import datetime
    import weathergen

    t0 = datetime(2022,6,1).timestamp()
    t1 = datetime(2022,9,1).timestamp()

    gen_times = np.arange(t0, t1, 60)

    weather = weathergen.generate(site='princeton', time=gen_times)

Note that the supplied year is arbitrary: the underlying model considers only annual and diurnal climatological variations. The supported sites are listed below, and are also stored in ``weathergen.sites``. Specified times should be supplied in Unix time.

The weather parameters are contained in the attributes of the ``weather`` object, e.g. ``weather.air_temp`` or ``weather.pressure``. The values are typically two-dimensional with shape ``(n_times, n_heights)``, where ``weather.time`` and ``weather.height`` describe the time and height of each dimension, in Unix time and meters above sea level. Single-level parameters are described by a single number for each time and do not have a layer dimension. 


Sites
-----

Supported sites are outlined below. Sites are included for the presence of major astronomical observatories, because they are climatological outliers with weird weather (and thus a good test of my model), or because they are cities that I like. 

.. list-table:: sites
   :widths: 25 50 50 75 30 30 30
   :header-rows: 1

   * - tag
     - description
     - country
     - notes
     - latitude (°N)
     - longitude (°E)
     - altitude (masl)
   * - tag
     - description
     - country
     - notes
     - latitude (°N)
     - longitude (°E)
     - altitude (masl)
