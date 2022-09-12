import matplotlib.pyplot as plt
from SWN_functions import *

def test_SWN():
    time_log = np.linspace(2, 6, 10**5)
    M_ej = 9.3 * u.solMass
    R = 624 * u.solRad
    E_exp = 10**51 * u.erg
    theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    plt.plot(time_log, L_obs(theta, (np.power(10, time_log)*u.s).to(u.d)))
    plt.grid()
    plt.yscale('log')
    plt.ylabel('L[erg/s]')
    plt.xlabel('t_log[log(s)]')
    plt.show()

if __name__ == '__main__':
    test_SWN()
    # M_ej = 9.3 * u.solMass
    # R = 624 * u.solRad
    # E_exp = 10 ** 51 * u.erg
    # theta = [R_x(R, 500).value, M_x(M_ej, 15).value, E_x(E_exp, 51).value, 0]
    # print(t_0(theta).to_value(u.s))
    # print(np.log10(180))