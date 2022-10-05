Weathergen
==========

``weathergen`` simulates a time-varying vertical profile of weather parameters. 

Usage
-----
       
After importing the package, simply pass a site and an array of times to the ```generate``` function. To simulate weather for the Chajnantor Plateau between June 1st and September 1st at a resolution of a minute, for example, we would write 

.. code-block:: python

    import numpy as np
    from datetime import datetime
    import weathergen

    t0 = datetime(2022,6,1).timestamp()
    t1 = datetime(2022,9,1).timestamp()

    gen_times = np.arange(t0, t1, 600)

    weather = weathergen.generate(site='chajnantor', time=gen_times)

Note that the supplied year is arbitrary; the underlying model considers only annual and diurnal climatological variations. The supported sites are listed below, and are also stored in ``weathergen.sites``. Specified times should be supplied in Unix time.

The weather parameters are contained in the attributes of the ``weather`` object, e.g. ``weather.air_temp`` or ``weather.pressure``. The values are typically two-dimensional with shape ``(n_times, n_heights)``, where ``weather.time`` and ``weather.height`` describe the time and height of each dimension, in Unix time and meters above sea level. Some parameters like ``weather.pwv`` and ``weather.cloud_cover`` are described by a single number for each time and do not have a layer dimension. 

Methodology
-----------

See paper. 

Sites
-----

Supported sites are shown below. Sites are chosen for the presence of astronomical observatories, or because I think that they're climatologically interesting.


