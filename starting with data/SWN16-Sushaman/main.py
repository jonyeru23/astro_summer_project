import matplotlib.pyplot as plt
from light_functions import *


M_ej = 9.3 * u.solMass
R = 624 * u.solRad
E_exp = 10**51 * u.erg
theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]

def test_mag_figs(band):
    time_log = np.linspace(2, 6)
    ob = S.ObsBandpass(band)
    plt.plot(time_log, [luminosity_to_mag(source_luminosity(theta, t, ob), 'ABMag') for t in (10**time_log * u.s).to(u.d)])
    plt.grid()
    plt.ylabel('M_AB')
    plt.xlabel('log time [log(s)]')
    plt.gca().invert_yaxis()
    plt.show()

def test_fig12():
    nus = np.linspace(14, 17)
    plt.plot(nus, [eq_2(theta, t_rc(theta), 10**nu) for nu in nus])
    plt.yscale('log')
    plt.grid()
    plt.xlabel('log(nu) [log(Hz)]')
    plt.ylabel('L_nu [erg/s/Hz]')
    plt.title('No 3 and no broken power law')
    plt.show()
    print(T_obs(theta, t_rc(theta)))


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


def test_fiter():
    v = S.ObsBandpass('v')
    plt.plot(v.wave, v.throughput)
    print(v.avgwave())
    print(v.equivwidth())
    print(v.rectwidth())
    plt.show()



if __name__ == '__main__':
    test_fig12()