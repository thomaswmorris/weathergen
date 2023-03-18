import numpy as np
from datetime import datetime

import weathergen

t0 = datetime(2023,1,1,1,0,0).timestamp()
t1 = datetime(2024,1,1,1,0,0).timestamp()

t_gen = np.arange(t0, t1, 3600)

def test_generation():

    weather = weathergen.Weather(region='princeton', verbose=True, generate=True)

    weather = weathergen.Weather(region='chajnantor', verbose=True, generate=False)
    weather.generate(time=t_gen, mode='random')

    weather = weathergen.Weather(region='tibet', verbose=True, generate=False)
    weather.generate(time=t_gen, mode='median')
