from SWN_functions import *
import pysynphot as S

h = const.h.to_value(u.erg / u.Hz)
k = const.k_B.to_value(u.erg / u.K)
distance = 26.4 * 10**6 * u.pc


def L_nu(theta, t, nu):
    temp_limit = (h * nu) / (3 * k * T_col(theta, t, nu))
    if temp_limit < 0.3:
        return Rayleign_Jeans(theta, t, nu)

    elif 0.3 < temp_limit < 3:
        return eq_2(theta, t, nu)


def T_col(theta, t, nu):
    return T_obs(theta, t) * (h * nu / (3 * k * T_obs(theta, t)))**0.2


def eq_2(theta, t, nu):
    return 0.9 * L_obs(theta, t) * 15 / np.pi**4 * \
           (h / (3 * k * T_col(theta, t, nu)))**4 * nu**3 * \
           (np.exp((h * nu)/(k * T_col(theta, t, nu))) - 1)**-1


def Rayleign_Jeans(theta, t, nu):
    return mag_to_flux(blackbody_mag(theta, t, nu))

def blackbody_mag(theta, t, nu):
    bb = S.BlackBody(T_col(theta, t, nu))

    bp = S.ObsBandpass(get_filter(nu))

    obs = S.Observation(bb, bp)

    mag = obs.effstim('VegaMag')

    return mag - 2.5 * np.log10(((R(theta, t) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / distance) ** 2))

def get_filter(nu):
    pass

def R(theta, t):
    """
        gets t in days, then converts the result to rsun
        R^2 = sqrt(L(t)/4pi sigma T^4)
        """
    sigma_sb = const.sigma_sb.to(u.erg * u.cm ** -2 * u.s ** -1 * u.K ** -4)
    return np.sqrt((L_obs(theta, t)*(u.erg/u.s) / (4 * np.pi * sigma_sb * (T_obs(theta, t)*u.k) ** 4))).to(u.solRad)


def mag_to_flux(mag):
    return 10 ** (-0.4 * mag)
