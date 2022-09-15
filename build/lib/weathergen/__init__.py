import numpy as np
from datetime import datetime

#f = resources.read_binary("weathergen", "generator.pkl")

import os
this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, 'generator.npy')


generator = np.load(DATA_PATH, allow_pickle=True)[()]

def get_tofd(t):
    DT = datetime.fromtimestamp(t)
    return DT.hour / 24 + DT.minute / 1440 + DT.second / 86400

def get_tofy(t):
    return datetime.fromtimestamp(t).timetuple().tm_yday / (366 + 1e-6)

def generate(region, t):
    
    i_tofy = np.digitize(get_tofy(t), bins=generator['tofy_bins']) - 1
    i_tofd = np.digitize(get_tofd(t), bins=generator['tofd_bins']) - 1

    gen   = generator[region]['gen'][i_tofy,i_tofd]
    data  = np.matmul(gen, np.random.standard_normal(gen.shape[1]))
    data *= generator[region]['rms'][i_tofy,i_tofd]
    data += generator[region]['avg'][i_tofy,i_tofd]

    weather = {}
    for col, d in zip(['lpwv', 'lbsd', 'lwvd', 'temperature', 'wind_east', 'wind_north'],
                            [data[0], data[1], *np.split(data[2:],4)]):
        weather[col] = d
        
    weather['pwv_rms']       = np.sqrt(np.exp(weather['lbsd']))
    weather['pwv']           = np.exp(weather['lpwv'])
    weather['water_density'] = np.exp(weather['lwvd'])
    weather['height']        = generator[region]['height']
    
    return weather