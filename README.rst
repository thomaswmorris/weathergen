Weathergen
==========

``weathergen`` generates time-varying weather parameters. 

.. list-table:: Single-level parameters
   :widths: 25 25 50 25
   :header-rows: 1

   * - Parameter
     - Tag
     - Description
     - Units
   * - Total column water vapor
     - total_water_vapor
     - The total amount of water vapor 
     - mm
   * - Total cloud cover
     - total_cloud_cover
     - The proportion of sky covered by clouds 
     - [none]
   * - Precipitation rate
     - total_precipitation
     - The rate of precipitation
     - mm/hr
     
.. list-table:: Multi-level parameters
   :widths: 25 25 50 25
   :header-rows: 1

   * - Parameter
     - Tag
     - Description
     - Units
   * - Air temperature 
     - air_temp
     - The temperature of the air 
     - K
   * - Air pressure
     - pressure
     - The air pressure
     - hPa

Usage
-----

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
    
Methodology
-----------

See paper. 

