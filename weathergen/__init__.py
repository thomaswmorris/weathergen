import numpy as np
import scipy as sp
import pandas as pd
import time, h5py, os, glob, re
from datetime import datetime
from scipy import interpolate

base, this_filename = os.path.split(__file__)

base = '/users/tom/desktop/repos/weathergen/weathergen'

def log_scale_invariant_spectrum(k,a,k0,nu):
    return np.log(a*(1+np.square(k/k0))**(-nu))

def get_sites(base):
    
    filenames = sorted(glob.glob(f'{base}/*.h5'))
    site_df = pd.DataFrame(columns=['name', 'longitude (°)', 'latitude (°)', 'altitude (m)'])
    for fn in filenames:
        with h5py.File(fn, 'r') as f:
            site_df.loc[len(site_df)] = re.findall(r'/([a-z_]+).h5', fn)[0], *[f[k][()] for k in ['lat', 'lon', 'alt']] 
    return site_df

sites = get_sites(base)
            
def get_day_hour(t):
    dt = datetime.fromtimestamp(t)
    return dt.hour + dt.minute/60 + dt.second/3600

def get_year_day(t):
    tt = datetime.fromtimestamp(t).timetuple()
    return (tt.tm_yday + tt.tm_hour / 24 - 1) 

def generate(t=time.time()+np.arange(0,86400,60), site='atacama'):
    
    if not site in sites.name.values:
        raise ValueError(f'\'{site}\' is not supported. Supported sites are:\n{sites}')
    
    filename = f'/users/tom/desktop/repos/weathergen/weathergen/{site}.h5'
    gen_data = {}
    with h5py.File(filename, 'r') as f:
        for key in list(f.keys()): gen_data[key] = f[key][()]

    if np.isscalar(t): raise ValueError('\'t\' must be not be a scalar (e.g. an array)')

    if np.array(t).ptp() < 30 * 86400 or not np.isclose(np.gradient(t), np.gradient(t).mean()).all():
        t_gen = np.mean(t) + np.arange(-15*86400, 15*86400, 60)
    else:
        t_gen = t.copy()
        
    n_gen, d_gen = len(t_gen), np.gradient(t_gen).mean()
    f_gen  = np.fft.fftfreq(n_gen, d_gen)
    yd_gen = list(map(get_year_day, t_gen))
    dh_gen = list(map(get_day_hour, t_gen))
    
    GLAFT = np.c_[[log_scale_invariant_spectrum(np.abs(f_gen), *p) for p in gen_data['laft_pars']]]
    
    GV = np.real(np.fft.ifft(4*np.exp(GLAFT)*np.sqrt(n_gen/d_gen)*np.exp(2j*np.pi*np.random.uniform(size=GLAFT.shape))))
    GD = np.matmul(gen_data['eigenmodes'][:,:GV.shape[0]], GV)

    _RGI = sp.interpolate.RegularGridInterpolator

    yd_points = gen_data['year_day_edge_points']
    dh_points = gen_data['day_hour_edge_points']
    binned_mean_data = gen_data['binned_mean_grid']
    binned_norm_data = gen_data['binned_norm_grid']

    GD *= _RGI((yd_points, dh_points), binned_norm_data, method='linear')((yd_gen, dh_gen)).T
    GD += _RGI((yd_points, dh_points), binned_mean_data, method='linear')((yd_gen, dh_gen)).T

    weather = {}
    for c, cg in zip(['lpwv', 'lwvd', 'temperature', 'wind_east', 'wind_north'], 
                     [GD[0], *np.split(GD[1:],4)]):

        weather[c] = cg
        
    weather['pwv']           = np.exp(weather.pop('lpwv'))
    weather['water_density'] = np.exp(weather.pop('lwvd'))
    weather['height']        = gen_data['hgt']

    for key in ['pwv', 'water_density', 'temperature', 'wind_east', 'wind_north']:

        weather[key] = sp.interpolate.interp1d(t_gen, weather[key])(t)
    
    return weather

