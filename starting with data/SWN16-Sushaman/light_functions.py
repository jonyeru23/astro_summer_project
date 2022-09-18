from SWN_functions import *
import pysynphot as S

h = const.h.to_value(u.erg / u.Hz)
k = const.k_B.to_value(u.erg / u.K)
distance = 26.4 * 10**6 * u.pc

"""
eventually I decided to leave eq2 as is, and make the wave length into freq
"""


def get_mag(theta, t, band_filter, system='VegaMag'):
    """
    the main function for retrieving the magnitude, devided to instences by the paper
    """
    ob = S.ObsBandpass(band_filter)
    nu = length_to_frequency(ob.avgwave())
    temp_limit = (h * nu) / (3 * k * T_col(theta, t, nu))
    if temp_limit < 0.3:
        return blackbody_mag(theta, t, ob, system)

    elif 0.3 < temp_limit < 3:
        return luminosity_to_mag(source_luminosity(theta, t, ob), system)

def get_range_length(ob, steps):
    """
    accepts the observation by the filter and the number of steps for the integral and outputs the relevant interval
    for the integral
    """
    start = ob.avgwave() - ob.equivwidth()/2
    stop = ob.avgwave() + ob.equivwidth()/2
    step_size = ob.equivwidth() / steps
    return np.arange(start, stop, step_size)


def get_range_freq(ob, steps):
    start = length_to_frequency(ob.avgwave() + ob.equivwidth()/2)
    stop = length_to_frequency(ob.avgwave() - ob.equivwidth()/2)
    step_size = (stop - start) / steps
    return np.arange(start, stop, step_size)


def source_luminosity(theta, t, ob):
    """
    makes the mag of the star by eq 2, by freqs
    """
    freqs = get_range_freq(ob, 10**3)
    Luminosity_to_hz = np.array([eq_2(theta, t, freq) for freq in freqs])
    return integral(freqs, Luminosity_to_hz)

# def L_ltt(theta, t):
#     """
#     time in days, integrates over the light time travel, analytically!
#     """
#     offset = theta[3]
#     t_offset = t - (offset * u.d)
#     if t < 0:
#         pass
#     if t < t_0(theta):
#         return integradet_L(t, L_smaller_t_0, theta) - integradet_L(t-t_rc(theta), L_smaller_t_0, theta)
#     elif t_0(theta) < t < t_s(theta):
#         return integradet_L(t, L_smaller_t_0, theta) - integradet_L(t-t_rc(theta), L_smaller_t_0, theta)
#
#
# def integrated_L(t, L_f)

def integral(x, y):
    return sum(y[i]*(x[i+1] - x[i]) for i in range(len(x) - 1))


def length_to_frequency(wave):
    """
    accepts wave in Angstrom and returns frequency in Hz
    """
    wave = (wave * u.angstrom).to_value(u.m)
    return const.c.value / wave


def L_nu(theta, t, band_filter, system='VegaMag'):
    """
    from the paper
    """
    ob = S.ObsBandpass(band_filter)
    nu = length_to_frequency(ob.avgwave())
    temp_limit = (h * nu) / (3 * k * T_col(theta, t, nu))
    if temp_limit < 0.3:
        return Rayleign_Jeans(theta, t, nu, system)

    elif 0.3 < temp_limit < 3:
        return source_luminosity(theta, t, ob)


def T_col(theta, t, nu):
    """
    eq. 3
    """
    return T_obs(theta, t) * (h * nu / (3 * k * T_obs(theta, t)))**0.2


def eq_2(theta, t, nu):
    """
    eq 2 so far
    """
    return 0.9 * L_7(theta, t) * 15 / np.pi**4 * \
           (h / (k * T_col(theta, t, nu)))**4 * nu**3 * \
           (np.exp((h * nu)/(k * T_col(theta, t, nu))) - 1)**-1


def Rayleign_Jeans(theta, t, nu, system):
    """
    the flux of a black body radiation
    """
    return mag_to_flux(blackbody_mag(theta, t, nu, system))

def blackbody_mag(theta, t, ob, system):
    """
    you get it right?, you need T_col because of the paper
    """
    bb = S.BlackBody(T_col(theta, t, nu=length_to_frequency(ob.avgwave())))


    obs = S.Observation(bb, ob)

    mag = obs.effstim(system)

    return mag - 2.5 * np.log10(((R(theta, t) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / distance) ** 2))


def R(theta, t):
    """
    gets t in days, then converts the result to rsun
    R^2 = sqrt(L(t)/4pi sigma T^4)
    """
    sigma_sb = const.sigma_sb.to(u.erg * u.cm ** -2 * u.s ** -1 * u.K ** -4)
    return np.sqrt((L_obs(theta, t)*(u.erg/u.s) / (4 * np.pi * sigma_sb * (T_obs(theta, t)*u.K) ** 4))).to(u.solRad)


def mag_to_flux(mag):
    """
    defenition of mag
    """
    return 10 ** (-0.4 * mag)

def luminosity_to_mag(L, system):
    """
    gets the absolute magnitude of a star, given L [erg/s], m_AB - m_Vega = 0.02
    returns the system in question
    """
    mag = -2.5 * np.log10(L/const.L_bol0.to_value(u.erg / u.s))
    if system == 'VegaMag':
        return mag - 0.02
    return mag

