from SWN_functions import *

def L_obs_new(theta, t):
    """
    time in days, returns dimensionless L
    """
    offset = theta[3]
    t_offset = t - offset*u.d
    L0 = L_smaller_t_0(theta)
    if t < 0:
        return L_smaller_0(theta, t_offset)

    elif t <= t_0(theta):
        return L0

    elif t_0(theta) <= t <= t_s(theta):
        return L0 * (t_offset/t_0(theta))**(-4/3)

    elif t_s(theta) <= t <= t_rec(theta):
        return L0 * (t_s(theta)/t_0(theta))**(-4/3) * (t_offset/t_s(theta))**(-0.35)
    return 0

def T_obs_new(theta, t):
    """
    time in days, returns dimensionless T
    """
    offset = theta[3]
    t_offset = t - (offset * u.d)
    T0 = 4.3 * 10 ** 5 * multi(theta, -0.17, -0.52, 0.35)
    if t < t_0(theta):
        return T0

    elif t_0(theta) <= t < t_s(theta):
        return T0 * (t_offset/t_0(theta))**(-0.45)

    elif t_s(theta) <= t < t_c(theta):
        return T0 * (t_s(theta)/t_0(theta))**(-0.45) * (t_offset/t_s(theta))**(-0.35)

    elif t_c(theta) <= t < t_rec(theta):
        return T0 * (t_s(theta)/t_0(theta))**(-0.45) * (t_c(theta)/t_s(theta))**(-0.35) * (t/t_c(theta))**(-0.6)

    return 0
