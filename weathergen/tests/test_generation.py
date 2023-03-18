import numpy as np
from datetime import datetime
import weathergen

t0 = datetime(2022,6,1).timestamp() # June 1, 2022
t1 = datetime(2022,9,1).timestamp() # August 1, 2022

gen_times = np.arange(t0, t1, 3600) # three months of hourly data

def test_generation():

    weather = weathergen.Weather(region='princeton', seasonal=True, diurnal=True)
    weather.generate(time=gen_times, mode='random')
