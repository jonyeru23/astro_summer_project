import matplotlib.pyplot as plt
from SWN_functions import *

def test_fig10():
    time_log = np.linspace(1.8, 6.2, 10**3)
    M_ej = 9.3 * u.solMass
    R = 624 * u.solRad
    E_exp = 10**51 * u.erg
    theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    time_log_real = time_log[10**time_log - t_rc(theta).value>0]
    plt.plot(time_log_real, T_max(theta, (10**time_log_real*u.s).to(u.d)))
    plt.grid()
    plt.yscale('log')
    plt.ylabel('L[erg/s]')
    plt.xlabel('t_log[log(s)]')
    plt.show()



def test_fig10_2():
    top = lambda theta, t: (L_smaller_t_0(theta)**-2 + L_smaller_t_s(theta, t)**-2)**-0.5
    bot = lambda theta, t: L_smaller_t_s(theta, t) + L_smaller_t_rec(theta, t)
    time_log = np.linspace(2, 6, 10 ** 3)
    M_ej = 9.3 * u.solMass
    R = 624 * u.solRad
    E_exp = 10 ** 51 * u.erg
    theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    plt.plot(time_log, L_f(theta, (np.power(10, time_log) * u.s).to(u.d), top))
    plt.plot(time_log, L_f(theta, (np.power(10, time_log) * u.s).to(u.d), bot), 'r')
    plt.grid()
    plt.yscale('log')
    plt.ylabel('L[erg/s]')
    plt.xlabel('t_log[log(s)]')
    plt.show()

def test_fig11():
    time = np.linspace(-850, 0, 10 ** 3)
    M_ej = 9.3 * u.solMass
    R = 624 * u.solRad
    E_exp = 10 ** 51 * u.erg
    theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    plt.plot(time, L_obs(theta, (time*u.s).to(u.d)))
    plt.grid()
    plt.yscale('log')
    plt.ylabel('L[erg/s]')
    plt.xlabel('t[s]')
    plt.xlim(-1000, 2000)
    plt.show()


if __name__ == '__main__':
    test_fig10()
    # test_fig10_2()
    # test_fig11()
    # M_ej = 9.3 * u.solMass
    # R = 624 * u.solRad
    # E_exp = 10 ** 51 * u.erg
    # theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    # print(t_0(theta).to_value(u.s))
    # print(t_s(theta).to_value(u.h))
    # print(np.log10(180))
    # print((10**2.873 * u.s).to(u.h))
