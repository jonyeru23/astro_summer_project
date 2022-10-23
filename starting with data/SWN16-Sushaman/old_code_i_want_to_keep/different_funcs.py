from SWN_functions import *
from light_functions import *


def L_ltt_7(theta, t):
    offset = theta[3]
    t_offset = t - offset * u.d
    L0 = L_smaller_t_0(theta)
    if t_offset < 0:
        return L_smaller_0(theta, t_offset)
    elif t_offset.to_value(u.h) < 0.20734:
        return (integrated_L(theta, t_offset, 0, L0) ** -2 +
                integrated_L(theta, t_offset, -4/3, L_smaller_t_s(theta, t_offset)) ** -2) ** -0.5
    else:
        return integrated_L(theta, t_offset, -4/3, L_smaller_t_s(theta, t_offset)) +\
               integrated_L(theta, t_offset, -0.35, L_smaller_t_rec(theta, t_offset))

def L_ltt_obs(theta, t):
    offset = theta[3]
    t_offset = t - offset * u.d
    L0 = L_smaller_t_0(theta)
    if t < 0:
        x = np.linspace(t_offset-t_rc(theta), t_offset, 100)
        return integral(x, [L_smaller_0(theta, time) for time in x]).value
    elif t <= t_0(theta):
        return integrated_L(theta, t_offset, 0, L0)

    elif t_0(theta) <= t <= t_s(theta):
        return integrated_L(theta, t_offset, -4 / 3, L_smaller_t_s(theta, t_offset))

    elif t_s(theta) <= t <= t_rec(theta):
        return integrated_L(theta, t_offset, -0.35, L_smaller_t_rec(theta, t_offset))
    return 0

def integrated_L(theta, t, a, L):
    return big_F_L(theta, t, a, L) - big_F_L(theta, t - t_rc(theta), a, L)


def big_F_L(theta, t, a, L):
    return 2 * L / t_rc(theta) * (t / (a + 1) * (1 - t/t_rc(theta)) + t**2 / (t_rc(theta)*(a+2)))


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


def L_7_new(theta, t):
    """
    the function to creat smooth power law
    """
    offset = theta[3]
    t_offset = t - offset
    if t_offset < 0:
        return L_smaller_0(theta, t_offset)
    elif t_offset.to_value(u.h) < 0.20734:
        return (L_smaller_t_0(theta) ** -2 + L_smaller_t_s(theta, t_offset) ** -2) ** -0.5
    else:
        return L_smaller_t_s(theta, t_offset) + L_smaller_t_rec(theta, t_offset)




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
