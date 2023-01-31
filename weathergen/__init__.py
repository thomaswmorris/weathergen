import time as ttime
import numpy as np
import scipy as sp
import pandas as pd
import h5py, os
from . import utils

RGI = sp.interpolate.RegularGridInterpolator

base, this_filename = os.path.split(__file__)
#base = '/users/tom/desktop/repos/weathergen/weathergen'

def get_sites(base):
    sites = pd.read_hdf(f'{base}/sites.h5').fillna('')
    sites = sites.loc[[os.path.exists(f'{base}/site_data/{tag}.h5') for tag in sites.index]]
    return sites

sites = get_sites(base)

site_string = sites.to_string(columns=['location', 'country', 'latitude', 'longitude', 'altitude'])
        
def generate(**kwargs): return Weather(**kwargs)

class Weather():

    def __init__(self, site=None, time=None, verbose=False, mode='random'):

        # mode is one of 'random', 'median'

        if site is None: site = 'princeton' # if no site is supplied, use 'princeton'
        if time is None: time = ttime.time() # if no time is supplied, use current time

        if not site in sites.index:
            raise ValueError(f"The site \'{site}\' is not supported. Available sites are:\n\n{site_string}")
    
        if not mode in ['random', 'median']:
            raise ValueError(f"mode must be one of ['random', 'median']")

        self.site, self.time, self.verbose = site, np.atleast_1d(time), verbose
        self.entry = sites.loc[self.site]
        self.lat, self.lon, self.alt = sites.loc[self.site, ['latitude', 'longitude', 'altitude']]

        self.quantiles = self.Quantiles(np.linspace(0,1,101))
        self.generate(mode=mode)

        # define some stuff
        if 'S' in self.entry.type:

            self.rain_precipitation_rate[self.rain_precipitation_rate < 1e-3] = 0 # we ignore precipitation less than 0.01 mm/hr
            self.snow_precipitation_rate[self.snow_precipitation_rate < 1e-3] = 0
            self.quantiles.rain_precipitation_rate[self.quantiles.rain_precipitation_rate < 1e-3] = 0
            self.quantiles.snow_precipitation_rate[self.quantiles.snow_precipitation_rate < 1e-3] = 0
            
        if 'P' in self.entry.type:

            self.cloud_cover = np.minimum(1, np.maximum(0, self.cloud_cover))

            for k in ['water_vapor', 'ozone', 'ice_water', 'snow_water', 'rain_water', 'liquid_water']:
                setattr(self, f'column_{k}', np.trapz(getattr(self, k), self.height, axis=0)[None])

        self.wind_bearing      = np.degrees(np.arctan2(self.wind_east, self.wind_north) + np.pi)
        self.wind_speed        = np.sqrt(np.square(self.wind_east) + np.square(self.wind_north))

        self.relative_humidity = np.minimum(100, np.maximum(1, utils.AH_to_RH(self.temperature, self.water_vapor)))
        self.dew_point         = utils.get_dew_point(self.temperature, self.relative_humidity)
        
    def generate(self, mode):

        filename = f'{base}/site_data/{self.site}.h5'
        self._gen_data = {}
        with h5py.File(filename, 'r') as f:
            for key in list(f.keys()): 
                self._gen_data[key] = f[key][()] if not key == 'gen_labels' else f[key][()].astype(str)

        dt_gen = np.gradient(self.time).min()
        gen_time = np.arange(self.time.min(), self.time.max() + dt_gen, dt_gen)

        n_quantiles = len(self._gen_data['q'])
        n_features  = len(self._gen_data['eigenmodes'])

        n_gen = len(gen_time)
        f_gen = np.fft.fftfreq(n_gen, dt_gen)

        self.year_day = list(map(utils.get_utc_year_day, gen_time))
        self.day_hour = list(map(utils.get_utc_day_hour, gen_time))

        yd_edge_index = self._gen_data['year_day_edge_index']
        dh_edge_index = self._gen_data['day_hour_edge_index']
        yd_knots = self._gen_data['year_day_edge_points']
        dh_knots = self._gen_data['day_hour_edge_points']

        self._gen_data['quantiles']  = self._gen_data['normalized_quantiles'].astype(float)
        self._gen_data['quantiles'] *= self._gen_data['quantiles_scaling']
        self._gen_data['quantiles'] += self._gen_data['quantiles_offset']

        if self.verbose: print('Generating time-ordered distributions . . .')

        qu = sp.interpolate.RegularGridInterpolator((yd_knots, dh_knots), 
        self._gen_data['quantiles'][yd_edge_index][:,dh_edge_index], 
        method='linear')((self.year_day, self.day_hour)).reshape(n_features * n_gen, n_quantiles, order='F')

        if mode == 'random':
            if self.verbose: print('Generating seeds . . .')
            GOOD_F_BINS = ~np.isnan(self._gen_data['azdft_binned'])
            AZDFT       = np.c_[[sp.interpolate.interp1d(self._gen_data['azdft_freq'][g], azdft[g], fill_value=0, bounds_error=False, kind='cubic')(f_gen) 
                                                        for azdft, g in zip(self._gen_data['azdft_binned'], GOOD_F_BINS)]]

            GEN_DFT = 2 * AZDFT * np.sqrt(n_gen / dt_gen) * np.exp(1j*np.random.uniform(low=0,high=2*np.pi,size=AZDFT.shape))
            GEN_V   = np.real(np.fft.ifft(GEN_DFT, axis=1))
            GD      = np.matmul(self._gen_data['eigenmodes'], GEN_V)

        else:
            GD = np.zeros((self._gen_data['eigenmodes'].shape[0], n_gen))
            
        gu = np.swapaxes(GD, 0, -1).reshape(n_features * n_gen, order='F')

        bin_data = sp.interpolate.interp1d(self._gen_data['z'], np.arange(n_quantiles), 
                                        bounds_error=False, fill_value='extrapolate')(gu)

        bin_data = np.minimum(np.maximum(bin_data, 0), n_quantiles - 2)

        DATA  = qu[np.arange(len(gu)), bin_data.astype(int)]
        DATA += (bin_data % 1) * (qu[np.arange(len(qu)), bin_data.astype(int) + 1] - DATA)
        DATA  = np.swapaxes(DATA.reshape(n_gen, n_features, order='F'), 0, 1)

        QUANTILES = np.swapaxes(qu.reshape(n_gen, n_features, n_quantiles, order='F'), 0, 1)

        feature_splits = np.r_[0, np.cumsum(self._gen_data['gen_widths'].astype(int))]

        for c, shaped, f0, f1 in zip(self._gen_data['gen_labels'], 
                                     self._gen_data['gen_shaped'], 
                                     feature_splits[:-1],
                                     feature_splits[1:],
                                     ):                                
        
            _resampled = sp.interpolate.interp1d(gen_time, DATA[f0:f1])(self.time)
            _quantiles = sp.interpolate.interp1d(self._gen_data['q'], QUANTILES[f0:f1], axis=-1, fill_value='extrapolate', kind='linear')(self.quantiles.q)

            if not shaped: 
                _resampled = _resampled.reshape(-1)
                _quantiles = _quantiles.reshape(-1, self.quantiles.n_q)
            setattr(self.quantiles, c, _quantiles.mean(axis=-2))
            setattr(self, c, _resampled)

        if 'P' in self.entry.type: setattr(self, 'height', self._gen_data['height'])
        
    class Quantiles():

        def __init__(self, q):

            self.q   = q
            self.n_q = len(q)