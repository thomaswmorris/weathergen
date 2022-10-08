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
       
In Python, import the package and pass a site and an array of times to the ```generate``` function. For example, to simulate weather for Princeton, New Jersey between June 1st and September 1st at a resolution of a minute we would write 

.. code-block:: python

    import numpy as np
    from datetime import datetime
    import weathergen

    t0 = datetime(2022,6,1).timestamp()
    t1 = datetime(2022,9,1).timestamp()

    gen_times = np.arange(t0, t1, 60)

    weather = weathergen.generate(site='princeton', time=gen_times)

Note that the supplied year is arbitrary; the underlying model considers only annual and diurnal climatological variations. The supported sites are listed below, and are also stored in ``weathergen.sites``. Specified times should be supplied in Unix time.

The weather parameters are contained in the attributes of the ``weather`` object, e.g. ``weather.air_temp`` or ``weather.pressure``. The values are typically two-dimensional with shape ``(n_times, n_heights)``, where ``weather.time`` and ``weather.height`` describe the time and height of each dimension, in Unix time and meters above sea level. Single-level parameters are described by a single number for each time and do not have a layer dimension. 

Methodology
-----------

See paper. 

Sites
-----

Supported sites are shown below. Sites are chosen for the presence of astronomical observatories, or because I think that they're climatologically interesting.

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
   * - arecibo
     - Arecibo, Puerto Rico
     - usa
     - Arecibo Observatory
     - 18.344
     - -66.753
     - 498
   * - armazones
     - Cerro Armazones
     - chile
     - ELT, VLT
     - -24.59
     - -70.192
     - 3046
   * - barrow
     - Barrow, Alaska
     - usa
     - GML
     - 71.323
     - -156.611
     - 11
   * - basrah
     - Basrah
     - iraq
     - Highest extreme temperatures
     - 30.526
     - 47.776
     - 5
   * - bure
     - Plateau de Bure
     - france
     - NOEMA
     - 44.634
     - 5.908
     - 2552
   * - cambridge
     - Cambridge, Massachusetts
     - usa
     - Harvard University, MIT
     - 42.374
     - -71.111
     - 8
   * - chajnantor
     - Cerro Chajnantor
     - chile
     - ACT, ALMA, APEX, ASTE, FYST, SO, TAO
     - -22.985
     - -67.741
     - 5040
   * - danakil
     - Danakil Desert
     - ethiopia
     - Highest average temperatures
     - 13.392
     - 40.821
     - -125
   * - effelsberg
     - Effelsberg
     - germany
     - ERT
     - 50.524
     - 6.883
     - 319
   * - falklands
     - Falkland Islands
     - uk
     - Mild Southern Ocean climate
     - -51.892
     - -59.221
     - 31
   * - graham
     - Mount Graham, Arizona
     - usa
     - LBT, VATT
     - 32.702
     - -109.89
     - 3178
   * - granada
     - Pico Veleta, Granada
     - spain
     - IRAM
     - 37.066
     - -3.393
     - 2850
   * - green_bank
     - Green Bank, West Virginia
     - usa
     - GBT
     - 38.43
     - -79.84
     - 807
   * - honolulu
     - Honolulu, Hawaii
     - usa
     - The nicest weather in the world
     - 21.382
     - -157.993
     - 8
   * - kerguelen
     - Kerguelen Islands
     - france
     - Extreme Southern Ocean climate
     - -49.349
     - 70.219
     - 10
   * - london
     - London
     - uk
     - The worst weather in the world
     - 51.477
     - 0.0
     - 12
   * - lucknow
     - Lucknow
     - india
     - Highest extreme PWV
     - 26.85
     - 80.95
     - 121
   * - malta
     - Malta
     - malta
     - Mediterranean climate
     - 35.881
     - 14.449
     - 90
   * - mauna_kea
     - Mauna Kea, Hawaii
     - usa
     - Mauna Kea Observatory
     - 19.823
     - -155.475
     - 4205
   * - mcmurdo
     - McMurdo Bay, Antarctica
     - antarctica
     - McMurdo Station
     - -77.846
     - 166.668
     - 10
   * - murchison
     - Murchison, Western Australia
     - australia
     - MRO, SKA
     - -26.703
     - 116.671
     - 395
   * - narrabri
     - Narrabri, New South Wales
     - australia
     - ATCA
     - -30.313
     - 149.55
     - 237
   * - ngari
     - Ngari, Tibet
     - china
     - AliCPT
     - 32.33
     - 80.03
     - 5176
   * - nobeyama
     - Nobeyama Observatory, Nagano
     - japan
     - 45m, NMA
     - 35.942
     - 138.476
     - 1350
   * - north_cape
     - Northern Cape
     - south africa
     - HERA, MeerKAT, SKA
     - -30.721
     - 21.411
     - 1075
   * - owens
     - Owens Valley, California
     - usa
     - OVRO
     - 37.232
     - -118.295
     - 1222
   * - pachon
     - Cerro Pachón, Chile
     - chile
     - LSST
     - -30.245
     - -70.749
     - 2663
   * - princeton
     - Princeton, New Jersey
     - usa
     - Princeton University
     - 40.344
     - -74.661
     - 58
   * - puna
     - Puna de Atacama
     - argentina
     - LLAMA
     - -24.192
     - -66.475
     - 4820
   * - quibdo
     - Quibdó, Colombia
     - colombia
     - Highest average PWV
     - 5.692
     - -76.658
     - 43
   * - samoa
     - American Samoa
     - usa
     - GML
     - -14.247
     - -170.564
     - 42
   * - singapore
     - Singapore
     - singapore
     - Very consistent climate
     - 1.354
     - 103.812
     - 15
   * - socorro
     - Socorro, New Mexico
     - usa
     - VLA
     - 34.1
     - -107.6
     - 2120
   * - south_pole
     - South Pole
     - antarctica
     - BICEP2, GML, Keck, SPT
     - -90.0
     - 0.0
     - 2835
   * - summit
     - Summit Camp, Greenland
     - denmark
     - GML, Summit Station
     - 72.579
     - -38.46
     - 3126
   * - teide
     - Mount Teide, Tenerife
     - spain
     - Teide Observatory
     - 28.3
     - -16.51
     - 2390
   * - washington
     - Mount Washington, New Hampshire
     - usa
     - Very erratic weather
     - 44.271
     - -71.303
     - 1917
   * - yakutsk
     - Yakutsk, Siberia
     - russia
     - Lowest extreme temperatures
     - 62.03
     - 129.73
     - 95
