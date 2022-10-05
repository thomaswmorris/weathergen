# Overview

```weathergen``` simulates a time-varying vertical profile of weather parameters. 

## Usage

After importing the package, simply pass a site and an array of times to the ```generate``` function:

```
import numpy as np
from datetime import datetime
import weathergen

t0 = datetime(2022,6,1).timestamp()
t1 = datetime(2022,9,1).timestamp()

gen_times = np.arange(t0, t1, 600)

weather = weathergen.generate(site='chajnantor', t=gen_times)
```

The supported sites are listed below, and are also stored in ```weathergen.sites```. Specified times should be supplied in Unix time.

## Methodology

See paper. 

## Sites

Supported sites are shown below. Sites are chosen for the presence of astronomical observatories, or because I think that they're climatologically interesting.


