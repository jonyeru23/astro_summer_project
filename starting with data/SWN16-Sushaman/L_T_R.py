from astropy import units as u
from astropy import constants as const
import numpy as np
import pysynphot as S

h = const.h.to_value(u.erg / u.Hz)
k = const.k_B.to_value(u.erg / u.K)
distance = 26.4 * 10 ** 6 * u.pc

"""
This file contains the following classes: Time, Luminosity, Temp, Filtered Luminosity.
concentrated into one place.

i need to deal with the t - offset.
i do it in one place and one place only, in FilteredL.get_filtered_L. after that all functions assume that t is after 
offset.
"""


class Time:
    @staticmethod
    def multi(theta, x, y, z) -> float:
        R500, M15, E51, offset = theta
        return (M15 ** x) * (R500 ** y) * (E51 ** z)

    def t_0(self, theta) -> float:
        """t_0 is converted to days"""
        return (155 * self.multi(theta, 0.23, 1.39, -0.81) * u.s).to(u.d)

    def t_s(self, theta) -> float:
        """t_s in converted to days"""
        return (3.6 * self.multi(theta, 0.44, 1.49, -0.56) * u.h).to(u.d)

    def t_c(self, theta) -> float:
        """t_c in days"""
        return 3.2 * self.multi(theta, 0.97, 2.02, -1.19) * u.d

    def t_rec(self, theta) -> float:
        """t_rec in days"""
        return 16.6 * self.multi(theta, 0.22, 0.76, -0.43) * u.d

    @staticmethod
    def t_rc(theta) -> float:
        """t_rc in days"""
        R500 = theta[0]
        return ((R500 * 500 * u.solRad).to(u.meter) / const.c).to(u.d)


class T(Time):
    def col(self, theta, t, nu) -> float:
        """nu is frequency"""
        m = 12
        return self.obs(theta, t) * (1 + (h * nu / (3 * k * self.obs(theta, t))) ** (-0.2 * m)) ** -(1 / m)
        # return self.obs(theta, t) * (h * nu / (3 * k * self.obs(theta, t))) ** 0.2

    def obs(self, theta, t):
        """eq 8 in Shusman"""
        if t < self.t_0(theta):
            return self.T_smaller_t_0(theta)

        elif self.t_0(theta) <= t < self.t_s(theta):
            return self.T_smaller_t_s(theta, t)

        elif self.t_s(theta) <= t < self.t_c(theta):
            return self.T_smaller_t_c(theta, t)

        elif self.t_c(theta) <= t < self.t_rec(theta):
            return self.T_smaller_t_rec(theta, t)

        return 1

    def T_smaller_t_0(self, theta):
        return 4.3 * 10 ** 5 * self.multi(theta, -0.17, -0.52, 0.35)

    def T_smaller_t_s(self, theta, t):
        return float(self.T_smaller_t_0(theta) * (t / self.t_0(theta)) ** (-0.45))

    def T_smaller_t_c(self, theta, t):
        return self.T_smaller_t_0(theta) * (self.t_s(theta) / self.t_0(theta)) ** (-0.45) * (t / self.t_s(theta)) ** (-0.35)

    def T_smaller_t_rec(self, theta, t):
        return self.T_smaller_t_0(theta) * (self.t_s(theta) / self.t_0(theta)) ** (-0.45) * (self.t_c(theta) /
                                                                                             self.t_s(theta)) ** (-0.35) * \
                   (t / self.t_c(theta)) ** (-0.6)


class Wave:
    def __init__(self, band_filter='V', system='VegaMag'):
        self.band_filter = band_filter
        self.system = system
        self.ob = S.ObsBandpass(band_filter)
        self.central_freq = self.length_to_frequency(self.ob.avgwave())

    @staticmethod
    def length_to_frequency(wave_length):
        """
        accepts wave in Angstrom and returns frequency in Hz
        """
        wave_length = (wave_length * u.angstrom).to_value(u.m)
        return const.c.value / wave_length

    @staticmethod
    def frequency_to_length(frequency):
        """
        accepts frequency in Hz and returns wave length in angstrom
        """
        return (const.c.value / frequency * u.m).to_value(u.angstrom)


    def get_range_freq(self):
        """
        return the frequency by the self.ob.wave object
        """
        return np.flip(self.length_to_frequency(self.ob.wave))


class L(Time):
    def __init__(self, steps=100, band_filter='V', system='VegaMag'):
        self.T = T()
        self.steps = steps
        self.wave = Wave(band_filter, system)
        self.system = system
        self.band_filter = band_filter

    def light_travel_time(self, theta, t, L_function) -> float:
        """eq 1 in Kozyreva paper"""
        dx = self.t_rc(theta) / self.steps
        return 2 / self.t_rc(theta) * \
               np.trapz([self.ltt_integrant(theta, t, t_tag, L_function) for t_tag in
                         np.linspace(t - self.t_rc(theta), t, self.steps)],
                        dx=dx)

    def ltt_integrant(self, theta, t, t_tag, L_function) -> float:
        """the function inside the integral"""
        return L_function(theta, t_tag) * (1 - (t - t_tag) / self.t_rc(theta))

    def eq_2(self, theta, t, nu):
        """eq 2 in Kozyreva paper, after light travel time, nu is freq"""
        return 0.9 * self.broken_power_law(theta, t) * 15 / np.pi ** 4 * \
               (h / (k * self.T.col(theta, t, nu))) ** 4 * nu ** 3 * \
               (np.exp((h * nu) / (k * self.T.col(theta, t, nu))) - 1) ** -1

    def broken_power_law(self, theta, t) -> float:
        """eq 7"""
        if t < 0:
            return self.smaller_0(theta, t)
        elif t.to_value(u.h) < 0.20734:
            return (self.smaller_t_0(theta) ** -2 + self.smaller_t_s(theta, t) ** -2) ** -0.5
        else:
            return self.smaller_t_s(theta, t) + self.smaller_t_rec(theta, t)

    def obs(self, theta, t) -> float:
        """eq 5"""
        if t <= 0:
            return self.smaller_0(theta, t)

        elif t <= self.t_0(theta):
            return self.smaller_t_0(theta)

        elif self.t_0(theta) <= t <= self.t_s(theta):
            return self.smaller_t_s(theta, t)

        elif self.t_s(theta) <= t <= self.t_rec(theta):
            return self.smaller_t_rec(theta, t)
        return 1

    def smaller_t_0(self, theta) -> float:
        return 1.8 * 10 ** 45 * self.multi(theta, -0.65, -0.11, 1.37)

    def smaller_0(self, theta, t) -> float:
        ratio = t / self.t_0(theta)
        return self.smaller_t_0(theta) * np.exp((-0.35 * ratio ** 2) + (0.15 * ratio))

    def smaller_t_s(self, theta, t) -> float:
        return self.smaller_t_0(theta) * (t / self.t_0(theta)) ** (-4 / 3)

    def smaller_t_rec(self, theta, t) -> float:
        return (self.smaller_t_0(theta) * (self.t_s(theta) / self.t_0(theta)) ** (-4 / 3) * (t / self.t_s(theta)) ** (
            -0.35)).to_value()



class Magnitude(L):
    """the default is to give back an absolute magnitude"""

    def get_filtered_abs_mag(self, theta, t):
        """the eq after eq 5 in Kozyreva"""
        offset = theta[3]
        t_offset = t - offset
        temp_limit = (h * self.wave.central_freq) / (3 * k * self.T.col(theta, t_offset, self.wave.central_freq))
        use_func = None
        if temp_limit < 0.3:
            use_func = self.Rayleign_Jean_pseudo_flux

        # theoretically, i need it to be smaller than 3 but i don't think we will get there
        # and i don't want to make discontinuity
        elif 0.3 < temp_limit:
            use_func = self.absolute_filtered_pseudo_flux
        return self.to_mag_from_pseudo_flux(self.light_travel_time(theta, t, use_func))

    def absolute_filtered_pseudo_flux(self, theta, t):
        return self.to_pseudo_flux_from_mag(self.absolute_mag_filtered(theta, t))

    def absolute_mag_filtered(self, theta, t):
        freq_range = self.wave.get_range_freq()
        flux = np.array([self.eq_2(theta, t, nu) for nu in freq_range]) / (
                    4 * np.pi * self.R(theta, t).to_value(u.cm) ** 2)
        wave_range = self.wave.frequency_to_length(freq_range)

        # flipped them so it would be from small to large
        sp = S.ArraySpectrum(np.flip(wave_range), np.flip(flux), fluxunits='fnu')
        return self.get_abs_mag_from_source(sp, theta, t)

    def Rayleign_Jean_pseudo_flux(self, theta, t):
        return self.to_pseudo_flux_from_mag(self.Rayleign_Jeans_Mag(theta, t))

    def Rayleign_Jeans_Mag(self, theta, t):
        """
        absolute magnitude of a blackbody, after convertion.
        without considering the distance of the supernova for generality.
        """
        bb = S.BlackBody(self.T.col(theta, t, nu=self.wave.central_freq))
        return self.get_abs_mag_from_source(bb, theta, t)

    def get_abs_mag_from_source(self, source, theta, t):
        obs = S.Observation(source, self.wave.ob)

        mag = obs.effstim(self.system)
        return mag - 2.5 * np.log10(((self.R(theta, t) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / (10.0 * u.pc)) ** 2))

    def R(self, theta, t):
        """
        gets t in days, then converts the result to rsun
        R = sqrt(L(t)/4pi sigma T^4)
        for this one i need the bolometric luminosity! so i will use broken power law with ltt
        """
        sigma_sb = const.sigma_sb.to(u.erg * u.cm ** -2 * u.s ** -1 * u.K ** -4)
        return np.sqrt(
            (self.light_travel_time(theta, t, self.broken_power_law) * (u.erg / u.s) / (4 * np.pi * sigma_sb *
                                                                                        (self.T.obs(theta,
                                                                                                    t) * u.K) ** 4))).to(
            u.solRad)

    @staticmethod
    def to_pseudo_flux_from_mag(mag):
        return 10 ** (-0.4 * mag)

    @staticmethod
    def to_mag_from_pseudo_flux(pseudo_flux):
        return -2.5 * np.log10(pseudo_flux)

    @staticmethod
    def convert_absolute_to_apparent(M, distance):
        return M + 5 * np.log10(distance.to_value()) - 5

    @staticmethod
    def convert_apparent_to_absolute(m, distance):
        return m + 5 - 5 * np.log10(distance.to_value())

    def absolute_mag_to_fluxerr(self, mag, mag_err):
        """
        a bit of derivation of errors, check with iair if this is the right way!!
        """
        return abs(self.to_pseudo_flux_from_mag(mag) * (-0.4) * np.log(10) * mag_err)



class FilteredL(L):
    """not exactly ccording to SOLID, but i need a few things in this calss to make it work properly"""

    def __init__(self, steps=100, band_filter='V', system='VegaMag'):
        super().__init__(steps)
        self.wave = Wave(band_filter, system)
        self.system = system
        self.band_filter = band_filter

    def get_filtered_L(self, theta, t):
        """the eq after eq 5 in Kozyreva"""
        offset = theta[3]
        t_offset = t - offset
        temp_limit = (h * self.wave.central_freq) / (3 * k * self.T.col(theta, t_offset, self.wave.central_freq))
        if temp_limit < 0.3:
            return self.Rayleign_Jeans(theta, t_offset)

        elif 0.3 < temp_limit < 3:
            return self.integrated_freq_L(theta, t_offset)
        return 1

    def integrated_freq_L(self, theta, t):
        freq_range = self.wave.get_range_freq()
        return np.trapz(y=[self.eq_2(theta, t, nu) for nu in freq_range],
                        x=freq_range)

    def Rayleign_Jeans(self, theta, t):
        """
        the flux of a black body radiation
        """
        return self.mag_to_flux(self.blackbody_absolute_mag(theta, t))

    def blackbody_absolute_mag(self, theta, t):
        """
        absolute magnitude of a blackbody, after convertion.
        without considering the distance of the supernova for generality.
        """
        bb = S.BlackBody(self.T.col(theta, t, nu=self.wave.central_freq))
        obs = S.Observation(bb, self.wave.ob)

        mag = obs.effstim(self.system)
        # new_mag = mag - 2.5 * np.log10(((self.R(theta, t) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / (10.0 * u.pc)) ** 2))
        return mag - 2.5 * np.log10(((self.R(theta, t) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / (10.0 * u.pc)) ** 2))

    @staticmethod
    def convert_absolute_to_apparent(M, distance):
        return M + 5 * np.log10(distance.to_value()) - 5

    def R(self, theta, t):
        """
        gets t in days, then converts the result to rsun
        R = sqrt(L(t)/4pi sigma T^4)
        """
        sigma_sb = const.sigma_sb.to(u.erg * u.cm ** -2 * u.s ** -1 * u.K ** -4)
        return np.sqrt((self.light_travel_time(theta, t) * (u.erg / u.s) / (4 * np.pi * sigma_sb *
                                                                            (self.T.obs(theta, t) * u.K) ** 4))).to(
            u.solRad)

    def mag_to_flux(self, mag):
        """
        defenition of mag
        """
        if self.system == 'VegaMag':
            mag += 0.02

        return const.L_bol0.to_value(u.erg / u.s) * 10 ** (-0.4 * mag)

    def mag_to_fluxerr(self, mag, mag_err):
        """
        a bit of derivation of errors, check with iair if this is the right way!!
        """
        return const.L_bol0.to_value(u.erg / u.s) * abs(self.mag_to_flux(mag) * (-0.4) * np.log(10) * mag_err)

    def luminosity_to_absolute_mag(self, L):
        """
        gets the absolute magnitude of a star, given L [erg/s], m_AB - m_Vega = 0.02
        returns the system in question
        """
        mag = -2.5 * np.log10(L / const.L_bol0.to_value(u.erg / u.s))
        if self.system == 'VegaMag':
            return mag - 0.02
        return mag

    def eq_2(self, theta, t, nu):
        """eq 2 in Kozyreva paper, after light travel time, nu is freq"""
        return 0.9 * self.light_travel_time(theta, t) * 15 / np.pi ** 4 * \
               (h / (k * self.T.col(theta, t, nu))) ** 4 * nu ** 3 * \
               (np.exp((h * nu) / (k * self.T.col(theta, t, nu))) - 1) ** -1
