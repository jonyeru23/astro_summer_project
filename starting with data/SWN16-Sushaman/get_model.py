from L_T_R import *
import emcee
"""
theta = [(R/500)[solRad], (M_ej/15)[solMass], (E_exp/(10**51))[erg], t_offset[dyas]]
i need to check with iair about the flux issue.
"""


class LogPosterior(FilteredL):
    @staticmethod
    def log_prior(theta):
        R500, M15, E51, toffset_d = theta

        if 0.6 < M15 < 2 and 0.1 < R500 < 4 and 0.1 < E51 < 5 and 0.70342593 < toffset_d < 0.73376157:
            return 0.0  # log(1)
        else:
            return -np.inf  # log(0)

    def log_likelihood(self, time, meas_mag, meas_mag_err, theta):
        """not finished, i need to sort out all the stuff regarding the apparent and absolute mag"""
        meas_flux = self.mag_to_flux(meas_mag)
        meas_flux_err = self.mag_to_fluxerr(meas_mag, meas_mag_err)

        expected_flux = np.array([self.get_filtered_L(theta, t) for t in time])

        chi2 = self.chi2(meas_flux, meas_flux_err, expected_flux)

        normalization = self.normalization(meas_flux)

        return -0.5 * np.sum(normalization + chi2)

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
        return lp + self.log_likelihood(theta, time, meas_mag, meas_mag_err)


class Sampler(LogPosterior):
    def __init__(self, steps=100, band_filter='V', system='VegaMag', nwalkers=8, ndim=4):
        super().__init__(steps, band_filter, system)
        self.nwalkers = nwalkers
        self.ndim = ndim

    def set_walkers(self):
        """
        initialize walkers -  theta = [R500, M15, E51, offset]
        the guess is from the other student
        """
        initial_guess = [1, 2, 2.5, 0.01]
        initial_spread = [0.2, 0.4, 0.5, 0.2]
        return np.random.randn(self.nwalkers, self.ndim) * initial_spread + initial_guess

    def write_sampler(self, file_name, x, y, yerr, steps):
        initial_pos = self.set_walkers()

        backend = emcee.backends.HDFBackend(file_name)
        backend.reset(self.nwalkers, self.ndim)

        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_posterior, backend=backend, args=(x, y, yerr))

        sampler.run_mcmc(initial_pos, steps, progress=True)



















