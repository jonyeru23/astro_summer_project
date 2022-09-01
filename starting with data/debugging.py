## remember to add extinction to their data!!!

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import astropy
import pysynphot as S
import emcee

k = 1
M_c = 1
distance = 26.4 * 10**6 * u.pc
# theta = [v85, R13, M_e, offset]

def f_p(n, M_e):
    """
    returns the ratio between the enitial mass and the core mass.
    """
    if n == 3 / 2:
        return (M_e / M_c) ** 0.5
    elif n == 3:
        return 0.08 * (M_e / M_c)


def L(n, t, theta):
    """
    t - days, M_e - solMass, v85 - vs * 10**-8.5 * cm * s**-1, R13 - R_e * 10**-13 * cm
    """
    v85, R13, M_e, offset = theta
    t_offset = t - offset
    if n == 3 / 2:
        nums = [1.88, -0.086, 1.67, 0.8]
    else:
        nums = [1.66, -0.175, 4.57, 0.73]
    A = nums[0] * 10 ** 42
    B = (v85 / (f_p(n, M_e) * (M_e + M_c) * k))
    C = (v85 ** 2 * (R13 / k))
    D = (nums[2] / (19.5 * np.sqrt(k * M_e * v85 ** -1)))

    return A * ((B * t_offset ** 2) ** nums[1]) * C * np.exp(-(D * t_offset) ** nums[3]) * u.erg / u.s


def T(n, t, theta):
    """
    the color temp is given by this function.
    where t is in days and M is in solar masses.
    """
    v85, R13, M_e, offset = theta
    t_offset = t - offset
    print(t, t_offset)
    if n == 3 / 2:
        nums = [2.05, 0.027]
    else:
        nums = [1.96, 0.016]
    T = nums[0] * 10 ** 4 * (((v85 * t_offset) ** 2) / (f_p(n, M_e) * (M_e + M_c) * k)) ** nums[1] * (
                R13 / k) ** 0.25 * t_offset ** -0.5
    return T * u.K


def R(n, t, theta):
    """
    gets t in days, then converts the result to rsun
    R^2 = sqrt(L(t)/4pi sigma T^4)
    """
    sigma_sb = const.sigma_sb.to(u.erg * u.cm ** -2 * u.s ** -1 * u.K ** -4)
    return np.sqrt(L(n, t, theta) / (4 * np.pi * sigma_sb * T(n, t, theta) ** 4)).to(u.solRad)


def get_mag(n, t, theta):
    mag = []
    for temp in T(n, t, theta).value:
        bb = S.BlackBody(temp)

        bp = S.ObsBandpass('v')

        obs = S.Observation(bb, bp)

        mag.append(obs.effstim('VegaMag'))
    return np.array(mag) - 2.5 * np.log10(((R(n, t, theta) / (1 * u.solRad)) ** 2) * ((1000.0 * u.pc / distance) ** 2))


def log_prior(theta):
    """
    i used the values Iair gave me, plus the times we expect the explosion to occur (by their paper)
    """
    v85, R13, M_e, offset = theta
    R_e = R13 * 1e13
    v_e = v85 * 10 ** 8.5
    if R_e > 1e10 and R_e < 1e14 and v_e > 100e5 and v_e < 100000e5 and M_e > 0.005 and M_e < 1 and \
            offset > 0.70342593 and offset < 0.73376157:
        return 0  # log(1)
    else:
        return -np.inf  # log(0)


def log_likelihood(n, theta, t, meas_mag, meas_mag_err):
    meas_flux = mag_to_flux(meas_mag)
    meas_flux_err = mag_to_fluxerr(meas_mag, meas_mag_err)

    expected_flux = mag_to_flux(get_mag(n, t, theta))

    chi2 = get_chi2(meas_flux, meas_flux_err, expected_flux)

    normalization = get_normalization(meas_flux)

    return -0.5 * np.sum(normalization + chi2)


def mag_to_flux(mag):
    return 10 ** (-0.4 * mag)


def mag_to_fluxerr(mag, mag_err):
    return abs(mag_to_flux(mag) * (-0.4) * np.log(10) * mag_err)


def get_chi2(meas, meas_err, expected):
    return ((meas - expected) / meas_err) ** 2


def get_normalization(meas):
    return np.log(2 * np.pi * meas ** 2)


def log_postirior(theta, n, t, meas_mag, meas_mag_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(n, theta, t, meas_mag, meas_mag_err)

def set_walkers():
    """
    # initialize walkers -  theta = [v85, R13, M_e, offset]
    """
    initial_guess = [1.7, 0.4, 0.02, 0.715]
    initial_spread = [0.2, 0.2, 0.002, 0.01]

    return np.random.randn(32, 4) * initial_spread + initial_guess



path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\excel files\combined_data.xlsx"

data = pd.read_excel(path, sheet_name='after extinction')

t = np.array(data.loc[:, 'JD - 2457651.0[day]'])
meas_mag = np.array(data.loc[:, 'V[mag]'])
meas_mag_err = np.array(data.loc[:, 'error_V[mag]'])




n = 3
initial_pos = set_walkers()
nwalkers, ndim = initial_pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_postirior,
                                args=(n, t, meas_mag, meas_mag_err))
sampler.run_mcmc(initial_pos, 5000, progress=True);