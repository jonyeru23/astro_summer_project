import unittest
from L_T_R import *
from astropy import units as u


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



class TestL(unittest.TestCase):
    def setUp(self) -> None:
        self.L = L()
        self.theta = [812/500, 8.09/15, 1.1, 0]
        self.t = 1

    def test_smaller_t_rec(self):
        self.assertEqual(type(self.L.smaller_t_rec(theta=self.theta, t=self.t)), np.float64)



if __name__ == '__main__':
    unittest.main()
