import numpy as np
import scipy as sp
import pandas as pd
import time, h5py, os, glob, re
from datetime import datetime
from scipy import interpolate

base, this_filename = os.path.split(__file__)
#base = '/users/tom/desktop/repos/weathergen/weathergen'

def log_scale_invariant_spectrum(k,a,k0,nu):
    return np.log(a*(1+np.square(k/k0))**(-nu))

def get_sites(base):

    sites = pd.read_csv(f'{base}/sites.csv', index_col=0)
    sites['has_data'] = [os.path.exists(f'{base}/site_data/{site_tag}.h5') for site_tag in sites.tag]
    return sites

sites = get_sites(base)
            
def get_day_hour(t):
    dt = datetime.fromtimestamp(t)
    return dt.hour + dt.minute/60 + dt.second/3600

def get_year_day(t):
    tt = datetime.fromtimestamp(t).timetuple()
    return (tt.tm_yday + tt.tm_hour / 24 - 1) 

def generate(t=[time.time()], site='chajnantor'):
    
    if not site in sites.tag.loc[sites.has_data].values:
        raise ValueError(f'\'{site}\' is not supported. Supported sites are:\n{sites.tag}')
    
    filename = f'/users/tom/desktop/repos/weathergen/weathergen/site_data/{site}.h5'
    gen_data = {}
    with h5py.File(filename, 'r') as f:
        for key in list(f.keys()): gen_data[key] = f[key][()].astype(float)

    if np.isscalar(t): raise ValueError('\'t\' must be not be a scalar (e.g. an array)')

    if np.array(t).ptp() < 60 * 86400 or not np.isclose(np.gradient(t), np.gradient(t).mean()).all():
        t_gen = np.mean(t) + np.arange(-30*86400, 30*86400, 3600)
    else:
        t_gen = t.copy()
        
    n_gen, d_gen = len(t_gen), np.gradient(t_gen).mean()
    f_gen  = np.fft.fftfreq(n_gen, d_gen)
    yd_gen = list(map(get_year_day, t_gen))
    dh_gen = list(map(get_day_hour, t_gen))
    
    #GLAFT = np.c_[[log_scale_invariant_spectrum(np.abs(f_gen), *p) for p in gen_data['laft_pars']]]
    GLAFT = np.c_[[sp.interpolate.interp1d(gen_data['freq_data'], _laft, bounds_error=False, fill_value=-np.inf)(f_gen) for _laft in gen_data['laft_data']]]
    
    GV  = np.real(np.fft.ifft(np.exp(GLAFT)*np.exp(2j*np.pi*np.random.uniform(size=GLAFT.shape))))
    GV /= GV.std(axis=1)[:,None]
    GD  = np.matmul(gen_data['eigenmodes'][:,:GV.shape[0]], GV)

    _RGI = sp.interpolate.RegularGridInterpolator

    yd_points = gen_data['year_day_edge_points']
    dh_points = gen_data['day_hour_edge_points']
    binned_mean_data = gen_data['binned_mean_grid']
    binned_norm_data = gen_data['binned_norm_grid']
    binned_min_data  = gen_data['binned_min_grid']
    binned_max_data  = gen_data['binned_max_grid']

    D = np.zeros(GD.shape)
    for i in range(D.shape[0]):
        D[i] = sp.interpolate.interp1d(gen_data['gq'], gen_data['QD'][i], bounds_error=False, fill_value='extrapolate')(GD[i])

    D *= _RGI((yd_points, dh_points), binned_norm_data, method='linear')((yd_gen, dh_gen)).T
    D += _RGI((yd_points, dh_points), binned_mean_data, method='linear')((yd_gen, dh_gen)).T
    #D  = np.minimum(gen_data['data_max'][:,None], _RGI((yd_points, dh_points), binned_max_data, method='linear')((yd_gen, dh_gen)).T)
    #D  = np.maximum(gen_data['data_min'][:,None], _RGI((yd_points, dh_points), binned_min_data, method='linear')((yd_gen, dh_gen)).T)

    weather = {}
    for c, cg in zip(['log_pwv', 'log_q', 'log_p', 'temperature', 'wind_east', 'wind_north'], 
                     [D[0], *np.split(D[1:],5)]):

        weather[c] = cg
        
    weather['pwv']           = np.exp(weather.pop('log_pwv'))
    weather['water_density'] = np.exp(weather.pop('log_q'))
    weather['pressure']      = np.exp(weather.pop('log_p'))
    weather['height']        = gen_data['hgt']

    for key in ['pwv', 'water_density', 'pressure', 'temperature', 'wind_east', 'wind_north']:

        weather[key] = sp.interpolate.interp1d(t_gen, weather[key], bounds_error=False, fill_value='extrapolate')(t)
        weather['pctls'] = gen_data['pctls'].T
    
    return weather

