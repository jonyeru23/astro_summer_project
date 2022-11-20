import unittest
from L_T_R import *
from astropy import units as u
import time
from get_model import *


theta = [812 / 500, 8.09 / 15, 1.1, 0]

class TestWave(unittest.TestCase):
    def setUp(self) -> None:
        self.V = Wave()
        self.U = Wave('U')
        self.R = Wave('R')
        self.all = [self.V, self.U, self.R]

    def test_length_to_freq(self):
        assert abs((self.V.central_freq / 1e10) - 5.4601e+14 / 1e10) < 1
        assert abs((self.U.central_freq / 1e10) - 8.3309e+14 / 1e10) < 1
        assert abs((self.R.central_freq / 1e10) - 4.6123e+14 / 1e10) < 1

    def test_range_freq(self):
        for wave in self.all:
            freqs = wave.get_range_freq()
            for i in range(9):
                assert freqs[i + 1] - freqs[i] > 0
            np.testing.assert_allclose(wave.ob.wave, np.flip(wave.frequency_to_length(freqs)), rtol=1e-02)
            # print(len(freqs))


class TestTime(unittest.TestCase):
    def setUp(self) -> None:
        self.theta = [812/500, 8.09/15, 1.1, 0]
        self.time = Time()

    def test_multi(self):
        self.assertEqual(type(self.time.multi(self.theta, 1, 1, 1)), float)

    def test_t0(self):
        self.assertEqual(self.time.t_0(self.theta).unit, u.d)

    def test_ts(self):
        self.assertEqual(self.time.t_s(self.theta).unit, u.d)

    def test_tc(self):
        self.assertEqual(self.time.t_c(self.theta).unit, u.d)

    def test_t_rec(self):
        self.assertEqual(self.time.t_rec(self.theta).unit, u.d)

    def test_t_rc(self):
        self.assertEqual(self.time.t_rc(self.theta).unit, u.d)


class TestT(unittest.TestCase):
    def setUp(self) -> None:
        self.T = T()
        self.theta = [812/500, 8.09/15, 1.1, 0]
        self.t = 1 * u.d

    def test_T_smaller_t_0(self):
        self.assertEqual(type(self.T.T_smaller_t_0(self.theta)), float)

    def test_T_smaller_t_s(self):
        self.assertEqual(type(self.T.T_smaller_t_s(self.theta, self.t)), float)

    def test_T_smaller_t_c(self):
        self.assertEqual(type(self.T.T_smaller_t_c(self.theta, self.t)), float)

    def test_T_smaller_t_rec(self):
        self.assertEqual(type(self.T.T_smaller_t_rec(self.theta, self.t)), float)

    def test_T_obs(self):
        times = [0 * u.d, self.T.t_0(self.theta), self.T.t_s(self.theta), self.T.t_c(self.theta),
                 self.T.t_rec(self.theta)]
        for time in times:
            self.assertEqual(type(self.T.obs(self.theta, time)), float)

    def test_T_col(self):
        times = [0 * u.d, self.T.t_0(self.theta), self.T.t_s(self.theta), self.T.t_c(self.theta),
                 self.T.t_rec(self.theta)]
        for time in times:
            self.assertEqual(type(self.T.col(self.theta, time, nu=5.4601e+14)), np.float64)


class TestL(unittest.TestCase):
    def setUp(self) -> None:
        self.L = L()
        self.theta = [812/500, 8.09/15, 1.1, 0]
        self.t = 1 * u.d
        self.times = [0*u.d, self.L.t_0(self.theta), self.L.t_s(self.theta), self.L.t_rec(self.theta)]
        self.nu = 5.4601e+14

    def test_smaller_t_rec(self):
        self.assertEqual(type(self.L.smaller_t_rec(theta=self.theta, t=self.t)), float)

    def test_smaller_t_s(self):
        self.assertEqual(type(self.L.smaller_t_s(theta=self.theta, t=self.t)), float)

    def test_smaller_0(self):
        self.assertEqual(type(self.L.smaller_0(theta=self.theta, t=self.t)), float)

    def test_smaller_t_0(self):
        self.assertEqual(type(self.L.smaller_t_0(theta=self.theta)), float)

    def test_obs(self):
        for time in self.times:
            self.assertEqual(type(self.L.obs(self.theta, time)), float)

    def test_eq2(self):
        for time in self.times:
            self.assertEqual(type(self.L.eq_2(self.theta, time, self.nu)), np.float64)

    def test_broken_power_law(self):
        for time in self.times:
            self.assertEqual(type(self.L.broken_power_law(self.theta, time)), float)


    """still need to test the rest of the functions"""
    # def test_integrant(self):
    #     for i, time in enumerate(self.times):
    #         if i == len(self.times) - 1:
    #             break
    #         for t_tag in np.linspace(time, self.times[i+1]):
    #             self.assertEqual(type(self.L.ltt_integrant(self.theta, time, t_tag, L.eq_2)), np.float64)


class TestMagnitude(unittest.TestCase):
    def setUp(self) -> None:
        self.mag = Magnitude()
        self.theta = [812 / 500, 8.09 / 15, 1.1, 0]
        self.times = [0 * u.d, self.mag.t_0(self.theta), self.mag.t_s(self.theta), self.mag.t_rec(self.theta)]
        self.t = 1 * u.d

    def test_absolute_mag_filtered(self):
        for t in self.times:
            use_func = self.mag.absolute_filtered_pseudo_flux
            t1 = time.perf_counter()
            mag1 = self.mag.to_mag_from_pseudo_flux(self.mag.light_travel_time(theta, t, use_func))
            t2 = time.perf_counter()
            print(f"Shusman {t2 - t1}")

            use_func = self.mag.Rayleign_Jean_pseudo_flux
            t1 = time.perf_counter()
            mag1 = self.mag.to_mag_from_pseudo_flux(self.mag.light_travel_time(theta, t, use_func))
            t2 = time.perf_counter()
            print(f"Rayleign_Jean {t2 - t1}")


    def test_get_mag_from_source(self):
        for t in self.times:
            bb = S.BlackBody(5e6)
            self.assertEqual(type(self.mag.get_abs_mag_from_source(bb, self.theta, t)), float)

    def test_R(self):
        for t in self.times:
            self.assertEqual(self.mag.R(self.theta, t).unit, u.solRad)


    def test_apparent_to_absolute(self):
        for absolute_mag in range(-15, 15):
            self.assertEqual(self.mag.convert_apparent_to_absolute(
                self.mag.convert_absolute_to_apparent(absolute_mag, distance),distance), absolute_mag)

    def test_pseudo_flux_to_mag(self):
        for absolute_mag in range(-15, 15):
            self.assertAlmostEqual(self.mag.to_mag_from_pseudo_flux(
                self.mag.to_pseudo_flux_from_mag(absolute_mag)), absolute_mag)

    def test_pseudo_fluxerr(self):
        for absolute_mag in range(-15, 15):
            for mag_error in np.linspace(0.01, 1):
                self.assertGreaterEqual(self.mag.absolute_mag_to_fluxerr(absolute_mag, mag_error), 0)

    # def test_R(self):
    #     for time in self.times:
    #         print(self.mag.R(self.theta, time))

class TestIntegrator(unittest.TestCase):
    def setUp(self) -> None:
        self.integrator = Integrator()
        self.func_example = lambda x: np.exp(-(x ** 2))
        self.mag = Magnitude()
        self.a = 0
        self.b = 1
        self.theta = [812 / 500, 8.09 / 15, 1.1, 0]
        self.error = 0.001

    def test_actual_func(self):
        t = (self.mag.t_0(theta) + self.mag.t_s(theta)) / 2
        # print(self.integrator.get_steps(self.mag.ltt_integrant, (t - self.mag.t_rc(theta)).to_value(), t.to_value(),
        #                                 theta, t, self.mag.absolute_filtered_pseudo_flux))
        for epsab in [1.49e-08, 1.49e-04, 1.49e-01]:
            t1 = time.perf_counter()
            l = [scipy.integrate.quad(self.mag.ltt_integrant, (t - self.mag.t_rc(theta)).to_value(), t.to_value(),
                                       args=(theta, t, self.mag.absolute_filtered_pseudo_flux), epsabs=epsab) for _ in range(10)]
            t2 = time.perf_counter()
            print(f"epabs = {epsab}: {t2-t1}")

class TestLogPosterior(unittest.TestCase):
    def setUp(self) -> None:
        self.theta = [812 / 500, 8.09 / 15, 1.1, 0]
        self.bad_theta = [812 / 500, 8.09 / 15, 1.1, 0]
        self.logPost = LogPosterior()
        path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\excel files\combined_data.xlsx"

        x, y, yerr = Sampler.get_data(path, sheet_name='after ext')
        self.x = x
        self.y = y
        self.yerr = yerr

    def test_prior(self):
        self.assertEqual(self.logPost.log_prior(self.theta), 0)
        self.assertEqual(self.logPost.log_prior(self.bad_theta), -np.inf)

    def test_chi2(self):
        x = np.linspace(0, 100)
        expected = x ** 2
        meas = expected * np.random.normal(len(x))
        error = np.array([0.1]*len(x))
        self.assertEqual(type(self.logPost.chi2(meas, error, expected)), np.ndarray)

    def test_log_likelihood(self):
        likihood = self.logPost.log_likelihood(self.x, self.y, self.yerr, self.theta)
        print(likihood)
        self.assertEqual(type(likihood), np.float64)

if __name__ == '__main__':
    unittest.main()
