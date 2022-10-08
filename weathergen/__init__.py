import numpy as np
import scipy as sp

import pandas as pd
import time as ctime
import h5py, os, glob, re, pytz
from datetime import datetime
from scipy import interpolate

RGI = sp.interpolate.RegularGridInterpolator

base, this_filename = os.path.split(__file__)
#base = '/users/tom/desktop/repos/weathergen/weathergen'

def scale_invariant_spectrum(k,a,k0,nu):
    return a*(1+np.square(k/k0))**(-nu)

def get_sites(base):

    sites = pd.read_csv(f'{base}/sites.csv', index_col=0)
    sites = sites.loc[[os.path.exists(f'{base}/site_data/{tag}.h5') for tag in sites.index]]
    return sites

sites = get_sites(base)
            
def get_utc_day_hour(t):
    dt = datetime.fromtimestamp(t, tz=pytz.utc).replace(tzinfo=pytz.utc)
    return dt.hour + dt.minute / 60 + dt.second / 3600

def get_utc_year_day(t):
    tt = datetime.fromtimestamp(t, tz=pytz.utc).replace(tzinfo=pytz.utc).timetuple()
    return (tt.tm_yday + get_utc_day_hour(t) / 24 - 1) 

class generate():

    def __init__(self, site=None, time=None):

        if site is None: site = 'princeton' # if no site is supplied, use 'princeton'
        if time is None: time = ctime.time() # if no time is supplied, use current time

        if not site in sites.index:
            raise ValueError(f'\'{site}\' is not supported. Supported sites are:\n{list(sites.index)}')
    
        self.site, self.time = site, np.atleast_1d(time)
        self.lat, self.lon, self.alt = sites.loc[sites.index == self.site, ['lat', 'lon', 'alt']].values[0]
        
        for k, v in self.generate().items(): setattr(self, k, v)

    def generate(self):

        filename = f'{base}/site_data/{self.site}.h5'
        gen_data = {}
        with h5py.File(filename, 'r') as f:
            for key in list(f.keys()): gen_data[key] = f[key][()].astype(float)
        self.gen_data = gen_data
        
        dt_gen = np.gradient(self.time).min()

        MINIMUM_SIM_DURATION = 30 * 86400 # simulate at least 60 days

        #t_min = np.minimum(self.time.mean() - 0.5 * MINIMUM_SIM_DURATION, self.time.min() - 86400)
        #t_max = np.maximum(self.time.mean() + 0.5 * MINIMUM_SIM_DURATION, self.time.max() + 86400)

        gen_time = np.arange(self.time.min(), self.time.max(), dt_gen)
            
        n_gen  = len(gen_time)
        f_gen  = np.fft.fftfreq(n_gen, dt_gen)

        self.yd_gen = list(map(get_utc_year_day, gen_time))
        self.dh_gen = list(map(get_utc_day_hour, gen_time))
        self.log_NDFT_pars = gen_data['log_NDFT_params']
        
        #GLAFT = np.c_[[sp.interpolate.interp1d(gen_data['freq_data'], _laft, bounds_error=False, fill_value=-np.inf)(f_gen) for _laft in gen_data['laft_data']]]
        
        AFT     = np.c_[[scale_invariant_spectrum(np.abs(f_gen), *p) for p in np.exp(self.log_NDFT_pars)]]
        GEN_DFT = AFT * np.sqrt(n_gen / dt_gen) * np.exp(1j*np.random.uniform(low=0,high=2*np.pi,size=AFT.shape))
        GV      = np.real(np.fft.ifft(GEN_DFT))
        GD      = np.matmul(gen_data['eigenmodes'][:,:GV.shape[0]], GV)

        yd_points = gen_data['year_day_edge_points']
        dh_points = gen_data['day_hour_edge_points']
        self.binned_mean_data = gen_data['binned_mean_grid']
        self.binned_norm_data = gen_data['binned_norm_grid']

        GD *= RGI((yd_points, dh_points), self.binned_norm_data, method='linear')((self.yd_gen, self.dh_gen)).T
        self.offset = RGI((yd_points, dh_points), self.binned_mean_data, method='linear')((self.yd_gen, self.dh_gen)).T
        GD += self.offset

        self.GD = GD.T
        self.QD = gen_data['QD']
        
        #binned_min_data  = gen_data['binned_min_grid']
        #binned_max_data  = gen_data['binned_max_grid']

        D = np.zeros(GD.shape)
        for i in range(D.shape[0]):
            D[i] = sp.interpolate.interp1d(gen_data['gq'], self.QD[i], bounds_error=False, fill_value='extrapolate')(GD[i])

        
        #D  = np.minimum(gen_data['data_max'][:,None], _RGI((yd_points, dh_points), binned_max_data, method='linear')((yd_gen, dh_gen)).T)
        #D  = np.maximum(gen_data['data_min'][:,None], _RGI((yd_points, dh_points), binned_min_data, method='linear')((yd_gen, dh_gen)).T)

        #D = np.maximum(gen_data['data_min'][:,None], np.minimum(gen_data['data_max'][:,None], D))

        cols = ['total_water_vapor', 'total_cloud_cover', 'total_precipitation', 
                'log_abs_hum', 'log_ozone', 'air_temp', 'wind_east', 'wind_north', 
                'pressure', 'log_cloud_cover', 'divergence', 'vorticity']

        weather = {}
        for c, cg in zip(cols, [D[0], D[1], D[2], *np.split(D[3:], 9)]):
            weather[c] = cg
        
        for key in cols:
            # print(key, weather[key].shape, gen_time.shape)
            weather[key] = sp.interpolate.interp1d(gen_time, weather[key], bounds_error=False, fill_value='extrapolate')(self.time).T
            if 'log_' in key: weather[key[4:]] = np.exp(weather.pop(key))

        weather['pctls']  = gen_data['pctls']
        weather['height'] = gen_data['height']

        weather['precipitation_rate'] = np.maximum(0, weather['total_precipitation'] - 1e-9 * weather['total_water_vapor'])

        # weather['ozone'] *= 1.27537 * (273.15 / 1e3) * weather['pressure'] / weather['air_temp']
        
        return weather

