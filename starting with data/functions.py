import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import astropy
import pysynphot as S
import emcee
from multiprocessing import Pool
import corner
import os
import shutil


k = 1
M_c = 1
distance = 26.4 * 10**6 * u.pc
labels = ["v85", "R13", "M_e", "offset"]
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
    if 1e10 < R_e < 1e14 and 100e5 < v_e < 100000e5 and 0.005 < M_e < 1 and \
            0.70342593 < offset < 0.73376157:
        return 0  # log(1)
    else:
        return -np.inf  # log(0)


def log_likelihood(n, theta, t, meas_mag, meas_mag_err):
    rel_t, rel_mag, rel_mag_err = get_rel_data(t, meas_mag, meas_mag_err, n, theta)

    meas_flux = mag_to_flux(rel_mag)
    meas_flux_err = mag_to_fluxerr(rel_mag, rel_mag_err)

    expected_flux = mag_to_flux(get_mag(n, rel_t, theta))

    chi2 = get_chi2(meas_flux, meas_flux_err, expected_flux)

    normalization = get_normalization(meas_flux)

    return -0.5 * np.sum(normalization + chi2)

def get_rel_data(t, mag, mag_err, n, theta):
    rel_t = t[t > t_bigger_than(n, theta)]
    rel_mag = mag[(len(t) - len(rel_t)):]
    rel_mag_err = mag_err[(len(t) - len(rel_t)):]
    return rel_t, rel_mag, rel_mag_err

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


def set_walkers(nwalkers):
    """
    # initialize walkers -  theta = [v85, R13, M_e, offset]
    """
    initial_guess = [1.7, 0.4, 0.02, 0.715]
    initial_spread = [0.2, 0.2, 0.002, 0.01]

    return np.random.randn(nwalkers, 4) * initial_spread + initial_guess


def t_bigger_than(n, theta):
    v85, R13, M_e, offset = theta
    return 0.2 * R13 * v85 ** -1 * max(0.5, (R13 ** 0.4 * (f_p(n, M_e) * k * (M_e + M_c) ** -0.2 * v85 ** -0.7)))


def t_smaller_than(n, theta):
    v85, R13, M_e, offset = theta
    return 7.4 * (R13 / k) ** 0.55


def get_data(file_path, sheet_name):
    data = pd.read_excel(file_path, sheet_name=sheet_name)

    t = np.array(data.loc[:, 'JD - 2457651.0[day]'])
    meas_mag = np.array(data.loc[:, 'V[mag]'])
    meas_mag_err = np.array(data.loc[:, 'error_V[mag]'])

    return t, meas_mag, meas_mag_err


def get_sampler(file_name, n, nwalkers, steps, x, y, yerr):
    initial_pos = set_walkers(nwalkers)
    nwalkers, ndim = initial_pos.shape

    backend = emcee.backends.HDFBackend(file_name)
    backend.reset(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_postirior, backend=backend, args=(n, x, y, yerr))

    sampler.run_mcmc(initial_pos, steps, progress=True)


def get_histograms(sampler, file_path, labels, n):
    plt.clf()
    fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(4):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    if 'binning' in file_path:
        plt.title(f'n={n} binning')
    else:
        plt.title(f'n={n} no binning')
    plt.savefig(file_path)


def get_lines_dots(flat_samples, x, y, yerr, n, file_path, start=0.01, end=3.5, close_up=False):
    plt.clf()
    inds = np.random.randint(len(flat_samples), size=100)
    time = np.linspace(start, end)
    for ind in inds:
        theta = flat_samples[ind]
        if close_up:
            time = np.linspace(t_bigger_than(n, theta), 0.02)
        v85, R13, M_e, offset = theta
        plt.plot(time, get_mag(n, time + offset, theta), "C1", alpha=0.1)

    plt.errorbar(x - offset, y, yerr=yerr, fmt=".k", capsize=0)

    plt.xlabel(f'JD - {2457651.0+offset}[day]')
    plt.ylabel('V[mag]')
    if 'binning' in file_path:
        plt.title(f'n={n} binning')
    else:
        plt.title(f'n={n} no binning')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.savefig(file_path)


def get_corner(flat_samples, file_path, labels, n):
    fig = corner.corner(flat_samples, labels=labels)
    fig.savefig(file_path)


def take_h5_files(folder_name):
    # os.mkdir(folder_name)

    ## move the h5 to the other folder
    for file in os.listdir():
        if 'h5' in file:
            shutil.move(file, f'/{folder_name}')


