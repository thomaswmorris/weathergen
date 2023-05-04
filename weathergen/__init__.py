# JMJ

import time as ttime
import numpy as np
import scipy as sp
import pandas as pd
import h5py, os
from . import utils

g = 9.80665

base, this_filename = os.path.split(__file__)
regions = pd.read_hdf(f'{base}/regions.h5').fillna('')
regions = regions.loc[[os.path.exists(f'{base}/region_data/{tag}.h5') for tag in regions.index]]
        
class InvalidRegionError(Exception):
    def __init__(self, invalid_region):
        region_string = regions.to_string(columns=['location', 'country', 'latitude', 'longitude', 'altitude'])
        super().__init__(f"The region \'{invalid_region}\' is not supported. Available regions are:\n\n{region_string}")

class Quantiles():

    def __init__(self, q):

        self.q = q

class Weather():

    def __init__(self, region="princeton", diurnal=True, seasonal=True, altitude=None, levels=None):

        if not region in regions.index:
            raise InvalidRegionError(region)

        self.region = regions.loc[region]
        self.diurnal, self.seasonal = diurnal, seasonal
        self.altitude = altitude if altitude is not None else self.region.altitude

        self.levels = levels if levels is not None else np.arange(self.altitude, 45000, 500)

        self.quantiles = Quantiles(q=np.linspace(0,1,101))
        filename = f"{base}/region_data/{self.region.name}.h5"

        with h5py.File(filename, "r") as f:
            for key in f.keys():

                if key in ["labels", "label_sources", "variable_labels", "single_variables", "profile_variables"]:
                    setattr(self, key, f[key][:].astype(str))
                else:
                    setattr(self, key, f[key][:].astype(f[key][:].dtype))

            self.eigenmodes = self.eigenmodes.astype(float)
            
            self.varcdf  = self.normalized_quantiles.astype(float)
            self.varcdf *= self.quantiles_scaling
            self.varcdf += self.quantiles_offset

            del self.normalized_quantiles
            del self.quantiles_scaling
            del self.quantiles_offset

        if not seasonal:
            self.varcdf = self.varcdf.mean(axis=0)[None]
            self.year_day_edge_points = [0,366]
            self.year_day_edge_index = [0,-1]

        if not diurnal:
            self.varcdf = self.varcdf.mean(axis=1)[:, None]
            self.day_hour_edge_points = [0,24]
            self.day_hour_edge_index = [0,-1]

        
    def generate(self, time=None, mode="random", fixed_quantiles={}, verbose=False):

        if time is None: 
            time = np.atleast_1d(ttime.time()) # if no time is supplied, use current time

        self.time = time

        dt_gen = 1e0 if self.time.size == 1 else np.gradient(self.time).min() 
        gen_time = np.arange(self.time.min() - dt_gen, self.time.max() + dt_gen, dt_gen)

        n_gen = gen_time.size
        f_gen = np.fft.fftfreq(n_gen, dt_gen)

        self.year_day = list(map(utils.get_utc_year_day, gen_time))
        self.day_hour = list(map(utils.get_utc_day_hour, gen_time))

        if verbose: print("Generating time-ordered distributions …")

        qu = sp.interpolate.RegularGridInterpolator((self.year_day_edge_points, self.day_hour_edge_points), 
                                                     self.varcdf[self.year_day_edge_index][:,self.day_hour_edge_index], 
                                                     method="linear")((self.year_day, self.day_hour)).reshape(len(self.eigenmodes) * n_gen, len(self.q), order="F")

        if len(fixed_quantiles) > 0:

            for label, val in fixed_quantiles.items():
                if not (0 < val < 1):
                    raise ValueError(f"All values in fixed_quantiles must be between 0 and 1 ({label}_quantile = {val}).")

            fixed_variables = np.isin(self.variable_labels, list(fixed_quantiles.keys()))
            eigenmode_pseudoinverse = np.linalg.pinv(self.eigenmodes[fixed_variables])
            desired_zscore = [np.sqrt(2)*sp.special.erfinv(2*fixed_quantiles[var]-1) for var in self.variable_labels[fixed_variables]]
            seed = np.matmul(eigenmode_pseudoinverse, desired_zscore)
            GD = np.repeat(np.matmul(self.eigenmodes, seed)[:, None], n_gen, axis=1)

        elif mode == "median":
            GD = np.zeros((self.eigenmodes.shape[0], n_gen))

        elif mode == "random":
            if verbose: print("mode: RANDOM")
            GOOD_F_BINS = ~np.isnan(self.azdft_binned)
            AZDFT       = np.c_[[sp.interpolate.interp1d(self.azdft_freq[g], azdft[g], fill_value=0, bounds_error=False, kind="cubic")(f_gen) 
                                                        for azdft, g in zip(self.azdft_binned, GOOD_F_BINS)]]
            GEN_DFT = 2 * AZDFT * np.sqrt(n_gen / dt_gen) * np.exp(1j*np.random.uniform(low=0,high=2*np.pi,size=AZDFT.shape))
            GEN_V   = np.real(np.fft.ifft(GEN_DFT, axis=1))
            GD      = np.matmul(self.eigenmodes, GEN_V)

        else: 
            raise ValueError(f'The specific method must be one of "random" or "median", or fixed_quantiles must be non-empty.')
            
        gu = np.swapaxes(GD, 0, -1).reshape(len(self.eigenmodes) * n_gen, order='F')

        bin_data = sp.interpolate.interp1d(self.z, np.arange(self.q.size), 
                                        bounds_error=False, fill_value='extrapolate')(gu)

        bin_data = np.minimum(np.maximum(bin_data, 0), self.q.size - 2)

        DATA  = qu[np.arange(len(gu)), bin_data.astype(int)]
        DATA += (bin_data % 1) * (qu[np.arange(len(qu)), bin_data.astype(int) + 1] - DATA)
        DATA  = np.swapaxes(DATA.reshape(n_gen, len(self.eigenmodes), order='F'), 0, 1)

        #QUANTILES = np.swapaxes(qu.reshape(n_gen, len(self.eigenmodes), len(self.q), order='F'), 0, 1)

        if 'S' in self.region.type:
            for label in self.single_variables:

                variable_mask = self.variable_labels == label
                resampled_data = sp.interpolate.interp1d(gen_time, DATA[variable_mask])(self.time)
                setattr(self, label, resampled_data.reshape(self.time.size))

        if 'P' in self.region.type:
            for label in self.profile_variables:

                variable_mask = self.variable_labels == label
                resampled_data = sp.interpolate.interp1d(gen_time, DATA[variable_mask])(self.time)
                setattr(self, label, resampled_data.reshape(-1, self.time.size))

            self.pressure = np.repeat(self.pressure_levels[:, None], self.time.size, axis=1)

            u_bin = np.minimum((self.geopotential[..., None] - g * self.levels < 0).sum(axis=0), self.geopotential.shape[0] - 2) - 1
            dummy_index = np.repeat(np.arange(self.time.size)[:, None], self.levels.size, axis=1)
            z_base = self.geopotential[u_bin, dummy_index]
            u_part = (g * self.levels - z_base) / (self.geopotential[u_bin + 1, dummy_index] - z_base)

            for attr in ['pressure', *self.profile_variables]:
                interp_data = (1 - u_part) * getattr(self, attr)[u_bin, dummy_index] + u_part * getattr(self, attr)[u_bin + 1, dummy_index]
                setattr(self, attr, interp_data.T)

        # define some stuff
        if 'S' in self.region.type:

            self.rain_precipitation_rate[self.rain_precipitation_rate < 1e-3] = 0 # we ignore precipitation less than 1 um/hr
            self.snow_precipitation_rate[self.snow_precipitation_rate < 1e-3] = 0

            self.surface_wind_bearing      = np.degrees(np.arctan2(self.surface_wind_east, self.surface_wind_north) + np.pi)
            self.surface_wind_speed        = np.sqrt(np.square(self.surface_wind_east) + np.square(self.surface_wind_north))
            self.surface_relative_humidity = utils.absolute_to_relative_humidity(self.surface_temperature, self.surface_water_vapor)
            self.surface_relative_humidity = np.minimum(100, np.maximum(1, self.surface_relative_humidity))
            self.surface_dew_point         = utils.get_dew_point(self.surface_temperature, self.surface_relative_humidity)
            
        if 'P' in self.region.type:

            self.cloud_cover = np.minimum(1, np.maximum(0, self.cloud_cover))

            for k in ['water_vapor', 'ozone', 'ice_water', 'snow_water', 'rain_water', 'liquid_water']:
                setattr(self, f'column_{k}', np.trapz(getattr(self, k), self.levels, axis=0))

            self.wind_bearing      = np.degrees(np.arctan2(self.wind_east, self.wind_north) + np.pi)
            self.wind_speed        = np.sqrt(np.square(self.wind_east) + np.square(self.wind_north))
            self.relative_humidity = utils.absolute_to_relative_humidity(self.temperature, self.water_vapor)
            self.relative_humidity = np.minimum(100, np.maximum(1, self.relative_humidity))
            self.dew_point         = utils.get_dew_point(self.temperature, self.relative_humidity)
