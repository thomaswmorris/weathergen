import numpy as npfrom datetime import datetimefrom importlib import resourcesimport pickle#f = resources.read_binary("weathergen", "generator.pkl")

import os
this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, 'weather_generator.npy')

#print(DATA_PATH)

generator = np.load(DATA_PATH,allow_pickle=True)[()]

#print(type(generator))#generator = pickle.loads(f)def get_tofd(t):    DT = datetime.fromtimestamp(t)    return DT.hour / 24 + DT.minute / 1440 + DT.second / 86400def get_tofy(t):    return datetime.fromtimestamp(t).timetuple().tm_yday / (366 + 1e-6)def generate(region,time,method='random'):    weather_cols = ['lpwv','lnwd','temperature', 'wind_east', 'wind_north']    i_tofy = np.digitize(get_tofy(time),bins=generator['tofy_bins']) - 1    i_tofd = np.digitize(get_tofd(time),bins=generator['tofd_bins']) - 1    weather = {}        if method == 'random':            u    = generator[region]['gen'][i_tofy,i_tofd]        gen  = np.sum(u*np.random.standard_normal(u.shape[1]),axis=1)         gen *= generator[region]['rms'][i_tofy,i_tofd]        gen += generator[region]['avg'][i_tofy,i_tofd]                for col,gen_data in zip(weather_cols,[gen[0],*np.split(gen[1:],4)]):            weather[col] = gen_data    if method == 'mean':        gen = generator[region]['avg'][i_tofy,i_tofd]                for col,gen_data in zip(weather_cols,[gen[0],*np.split(gen[1:],4)]):            weather[col] = gen_data                weather['pwv'] = np.exp(weather['lpwv'])    weather['water_density'] = np.exp(weather['lnwd']+weather['lpwv'])    weather['height'] = generator[region]['height']        return weather