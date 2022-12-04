from L_T_R import *
import emcee
import pandas as pd
from multiprocessing import Pool
"""
theta = [(R/500)[solRad], (M_ej/15)[solMass], (E_exp/(10**51))[erg], t_offset[dyas]]
"""


class LogPosterior(Magnitude):
    @staticmethod
    def log_prior(theta):
        """using what iair gave me, for mass and offset, the rest from the other student"""
        R500, M15, E51, toffset_d = theta

        if 0.0005 < R500 < 4 and 0.6 < M15 < 2 and 0.1 < E51 < 5 and 0.70342593 < toffset_d < 0.73376157:
            return 0.0  # log(1)
        else:
            return -np.inf  # log(0)

    def log_likelihood(self, time, meas_mag, meas_mag_err, theta):
        """
        my data is apparent magnitude. i converted it to absolute magnitude and then to pseudo flux
        the function i am comparing it to gives me absolute magnitude, and i convert it to pseudo flux as well
        then i can compare the two and get chi2 and everything.

        all the data is in the form of np.array.
        """
        absolute_mag = self.convert_apparent_to_absolute(meas_mag, distance)
        absolute_mag_err = meas_mag_err

        meas_flux = self.to_pseudo_flux_from_mag(absolute_mag)
        meas_flux_err = self.absolute_mag_to_fluxerr(absolute_mag, absolute_mag_err)

        expected_flux = np.array([self.to_pseudo_flux_from_mag(self.get_filtered_abs_mag(theta, t)) for t in time])

        meas_flux, meas_flux_err, expected_flux = self.delete_beyond_validity(meas_flux, meas_flux_err, expected_flux)

        chi2 = self.chi2(meas_flux, meas_flux_err, expected_flux)

        normalization = self.normalization(meas_flux)

        return -0.5 * np.sum(normalization + chi2)

    @staticmethod
    def delete_beyond_validity(meas_flux, meas_flux_err, expected_flux):
        expected_flux = expected_flux[expected_flux > 1]
        meas_flux = meas_flux[:len(expected_flux)]
        meas_flux_err = meas_flux_err[:len(expected_flux)]
        return meas_flux, meas_flux_err, expected_flux

    @staticmethod
    def chi2(meas, meas_err, expected):
        """might take a single value or a np.array"""
        return ((meas - expected) / meas_err) ** 2

    @staticmethod
    def normalization(meas):
        """might take a single value or a np.array"""
        return np.log(2 * np.pi * meas ** 2)

    def log_posterior(self, theta, time, meas_mag, meas_mag_err):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(time, meas_mag, meas_mag_err, theta)


class Sampler(LogPosterior):
    def __init__(self, data_file_path, sheet_name, steps=100, band_filter='V', system='VegaMag', nwalkers=8, ndim=4):
        super().__init__(steps, band_filter, system)
        self.nwalkers = nwalkers
        self.ndim = ndim

        x, y, yerr = self.get_data(data_file_path, sheet_name)
        self.x = x
        self.y = y
        self.yerr = yerr


    @staticmethod
    def get_data(file_path, sheet_name):
        data = pd.read_excel(file_path, sheet_name=sheet_name)

        t = np.array(data.loc[:, 'JD - 2457651.0[day]'])
        meas_mag = np.array(data.loc[:, 'V[mag]'])
        meas_mag_err = np.array(data.loc[:, 'error_V[mag]'])

        return t, meas_mag, meas_mag_err

    def set_walkers(self):
        """
        initialize walkers -  theta = [R500, M15, E51, offset]
        the guess is from the other student
        """
        initial_guess = [1, 1, 2.5, 0.715]
        initial_spread = [0.2, 0.2, 0.5, 0.01]
        return np.random.randn(self.nwalkers, self.ndim) * initial_spread + initial_guess

    def write_sampler(self, sampler_file_name, steps):
        """
        the backend allows to write it to a file, for later use.
        important - the first argument to self.log_posterior should be theta so that the walkers cna walk properly
        """
        initial_pos = self.set_walkers()

        backend = emcee.backends.HDFBackend(sampler_file_name)
        backend.reset(self.nwalkers, self.ndim)
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_posterior, backend=backend,
                                            args=(self.x, self.y, self.yerr), pool=pool)

            sampler.run_mcmc(initial_pos, steps, progress=True)




















