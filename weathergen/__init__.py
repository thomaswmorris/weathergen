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

def generate(site,t):
    
    i_tofy = np.digitize(get_tofy(t),bins=generator['tofy_bins']) - 1
    i_tofd = np.digitize(get_tofd(t),bins=generator['tofd_bins']) - 1

    u    = generator[site]['gen'][i_tofy,i_tofd]
    gen  = np.sum(u*np.random.standard_normal(u.shape[1]),axis=1) 
    gen *= generator[site]['rms'][i_tofy,i_tofd]
    gen += generator[site]['avg'][i_tofy,i_tofd]

    weather = {}
    for col,gen_data in zip(['lpwv', 'lbsd', 'lwvd', 'temperature', 'wind_east', 'wind_north'],
                            [gen[0], gen[1], *np.split(gen[2:],4)]):
        weather[col] = gen_data
        
    weather['pwv_rms']       = np.sqrt(np.exp(weather['lbsd']))
    weather['pwv']           = np.exp(weather['lpwv'])
    weather['water_density'] = np.exp(weather['lwvd'])
    weather['height']        = generator[site]['height']
    
    return weather