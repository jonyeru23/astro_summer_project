from astropy import units as u
import numpy as np
"""
theta = [R500, M15, E51, offset]

All the functions will be designed to accept arrays, not individuals points
The default value of time is days.
R, M, E are dimensionless
Only time has units 
"""
###########################################################################
def L_obs(theta, time):
    """
    time in days
    """
    offset = theta[3]
    time_offset = time - offset
    L = np.zeros(time_offset.shape)
    for i, t in enumerate(time_offset):
        if t <= 0:
            ratio = t/t_0(theta)
            L[i] = 1.8 * 10**45 * multi(theta, -0.65, -0.11, 1.37) \
                   * np.exp(-0.35 * ratio**2 + 0.15*ratio)

        elif t <= t_0(theta):
            L[i] = 1.8 * 10**45 * multi(theta, -0.65, -0.11, 1.37)

        elif t_0(theta) <= t <= t_s(theta):
            L[i] = 2.7 * 10**43 * multi(theta, -0.34, 1.74, 0.29) \
                     * t.to_value(u.h)**(-4/3)

        elif t_s(theta) <= t <= t_rec(theta):
            L[i] = 1.6 * 10**42 * multi(theta, -0.78, 0.28, 0.84) \
                     * t.value**-0.35
    return L


def T_obs(theta, time):
    """
    time in days
    """
    offset = theta[3]
    time_offset = time - offset
    T = np.zeros(time_offset.shape)
    for i, t in enumerate(time_offset):
        if t < t_0(theta):
            T[i] = 4.3 * 10**5 * multi(theta, -0.17, -0.52, 0.35)

        elif t_0(theta) <= t < t_s(theta):
            T[i] = 10**5 * multi(theta, -0.07, 0.1, -0.01) \
                   * t.to_value(u.h)**-0.45

        elif t_s(theta) <= t < t_c(theta):
            T[i] = 3 * 10**4 * multi(theta, -0.11, -0.04, 0.04) * t.value**-0.35

        elif t_c(theta) <= t < t_rec(theta):
            T[i] = 4.1 * 10**4 * multi(theta, 0.13, 0.46, -0.25) * t.value**-0.6
    return T


##############################################################################
def t_0(theta):
    """t_0 is converted to days"""
    return (155 * multi(theta, 0.23, 1.39, -0.81) * u.s).to(u.d)


def t_s(theta):
    """t_s in converted to days"""
    return (3.6 * multi(theta, 0.44, 1.49, -0.56) * u.h).to(u.d)


def t_c(theta):
    """t_c in days"""
    return 3.2 * multi(theta, 0.97, 2.02, -1.19) * u.d


def t_rec(theta):
    """t_rec in days"""
    return 16.6 * multi(theta, 0.22, 0.76, -0.43) * u.d


def multi(theta, x, y, z):
    R500, M15, E51, offset = theta
    return (M15**x) * (R500**y) * (E51**z)


#########################################################################

def L_to_rel_mag(L):
    return -2.5*np.log10(L)


def M_x(M, x):
    return M/x


def R_x(R, x):
    return R/x


def E_x(E, x):
    return E/(10**x)










