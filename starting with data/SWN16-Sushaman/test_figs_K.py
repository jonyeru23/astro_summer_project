import matplotlib.pyplot as plt
from L_T_R import *

"""this file is for testing the figures from Kozyreva and a bit from Shusman"""

M_ej = 9.3 * u.solMass
R = 624 * u.solRad
E_exp = 10**51 * u.erg
theta = [(R/500).value, (M_ej/15).value, (E_exp/(10**51)).value, 0]

m12l15rot2 = [812/500, 8.09/15, 1.1, 0]
m12l5rot2 = [345/500, 9.15/15, 1.1, 0]
m15l15rot0 = [1024/500, 11.33/15, 1.1, 0]
m15l5rot8 = [268/500, 6.1/15, 1.1, 0]

progenitors = {
    'm12l15rot2': m12l15rot2,
    'm12l5rot2': m12l5rot2,
    'm15l15rot0': m15l15rot0,
    'm15l5rot8': m15l5rot8
}

def test_classes():
    Lumi = L(100)
    data = pd.read_csv(
        r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\SWN16-Sushaman\test_data\noraml_eta.csv")
    # print(data.loc[:, 'L_bol[erg/s]'])
    # print([L_obs(theta, (t*u.s).to(u.d)) for t in data.loc[:, 't[s]']])
    # print(type(L_obs_new(theta, (data.loc[0, 't[s]']*u.s).to(u.d)).value))
    plt.plot(np.log10(data.loc[:, 't[s]']), np.array(data.loc[:, 'L_bol[erg/s]']), 'r-', label='Test func')
    # plt.plot(np.log10(data.loc[:, 't[s]']), np.array([float(L_obs_new(theta, (t*u.s).to(u.d))) for t in data.loc[:, 't[s]']]), label='eq 5')
    plt.plot(np.log10(data.loc[:, 't[s]']),
             np.array([float(Lumi.broken_power_law(theta, (t * u.s).to(u.d))) for t in data.loc[:, 't[s]']]), label='broken law eq 7')
    plt.legend()
    plt.grid()
    # plt.xlim(-1, 100)
    plt.ylabel('L[erg/s]')
    plt.xlabel('log(t[s])')
    plt.title('R=624 solRad, M=9.3 solMass, E=1e51 erg')
    plt.show()


def test_fig_3():
    L_ltt = L(100).light_travel_time
    figure, axis = plt.subplots(2, 2)
    prog = list(progenitors.items())
    index = 0
    time = np.linspace(1, 6)
    for row in range(2):
        for cul in range(2):
            axis[row, cul].plot(time, [np.log10(L_ltt(prog[index][1], (t * u.s).to(u.d))) for t in 10**time])
            axis[row, cul].set_title(prog[index][0])
            axis[row, cul].grid()
            axis[row, cul].set_ylim(42, 46)
            index += 1
    plt.show()

def test_fig_5():
    log_nu = np.linspace(13, 18)
    L_nu = FilteredL(100).eq_2
    plt.plot(log_nu, [np.log10(L_nu(m12l15rot2, (300*u.s).to(u.d), nu)) for nu in 10**log_nu])
    plt.grid()
    plt.ylim(25, 29)
    plt.xlabel('log nu [log(Hz)]')
    plt.ylabel('log L_nu [log(erg / (sHz))')
    plt.title('Figure 5 Kozyreva')
    plt.show()

def test_fig_9():
    """there is a problem somewhere"""
    log_t = np.linspace(1.5, 6.5)
    to_mag = FilteredL(100).luminosity_to_absolute_mag
    L_filtered = FilteredL(100).get_filtered_L
    plt.plot(log_t, [to_mag(L_filtered(m12l15rot2, (t*u.s).to(u.d))) for t in 10**log_t])
    plt.show()


def test_fig_10_shusman():
    """here it works"""
    log_time = np.linspace(1.8, 6)
    Lumi = FilteredL(100)
    T_obs = Lumi.T.obs
    nu = Lumi.wave.central_freq
    plt.plot(log_time, [np.log10(T_obs(m12l15rot2, (t*u.s).to(u.d))) for t in 10**log_time])
    plt.grid()
    plt.xlabel('log(t[s])')
    plt.ylabel('log(T[K])')
    plt.title('Figure 10 Shusman')
    plt.show()

def test_fig_12():
    """so far no luck here, i dont know what is the problem"""
    log_time = np.linspace(0.5, 2)
    Lumi = FilteredL(100, 'U')
    T_col = Lumi.T.col
    nu = Lumi.wave.central_freq
    print(T_col(m12l15rot2, 1*u.d, nu))
    plt.plot(log_time, [T_col(m12l15rot2, t*u.d, nu) for t in 10**log_time])
    plt.grid()
    plt.xlabel('log(t[d])')
    plt.ylabel('T[K])')
    plt.title('Figure 12 Kyzoreva')
    plt.show()

def test_fig_4():
    """wwent well"""
    T_col = L().T.col
    nu = FilteredL(band_filter='V').wave.central_freq
    figure, axis = plt.subplots(2, 2)
    prog = list(progenitors.items())
    index = 0
    time = np.linspace(1, 6)
    for row in range(2):
        for cul in range(2):
            axis[row, cul].plot(time, [np.log10(T_col(prog[index][1], (t * u.s).to(u.d), nu)) for t in 10 ** time])
            axis[row, cul].set_title(prog[index][0])
            axis[row, cul].grid()
            axis[row, cul].set_ylim(4, 6)
            index += 1
    plt.show()

def test_fig_17_19_shusman(length):
    """
    length by the paper
    test passed,
    not exactly the same graph but i think it's because of the difference between eq 39 and 5.
    """
    log_time = np.linspace(1.8, 6.2)
    lumi = FilteredL(band_filter='R', system='ABMag')
    to_mag = lumi.luminosity_to_absolute_mag
    eq_2 = lumi.eq_2
    nu = lumi.wave.length_to_frequency((length * 10**-9 * u.m).to_value(u.angstrom))
    print(nu)
    # plt.plot(log_time, [to_mag(lumi.get_filtered_L(theta, t*u.s)) for t in 10**log_time])
    # plt.plot(log_time, [to_mag(eq_2(theta, (t*u.s).to(u.d), nu)*nu) for t in 10**log_time])
    plt.plot(log_time, [to_mag(lumi.integrated_freq_L(theta, (t*u.s).to(u.d))) for t in 10**log_time])
    plt.ylim(-12, -20)
    plt.xlabel('log t [s]')
    plt.ylabel('M_AB')
    plt.title(f'Lambda = {length}nm')
    plt.grid()
    plt.show()

def test_magnitude_graph():
    """my own test to see if the integral over the wave lengths is ok"""
    log_time = np.linspace(1.8, 6)
    lumi = FilteredL()
    to_mag = lumi.luminosity_to_absolute_mag
    to_apparent = lumi.convert_absolute_to_apparent
    plt.plot(log_time, [to_apparent(to_mag(lumi.get_filtered_L(theta, (t * u.s).to(u.d))), distance) for t in 10**log_time], 'k', label='filtered L')
    # plt.plot(log_time, [to_mag(lumi.integrated_freq_L(theta, (t * u.s).to(u.d))) for t in 10 ** log_time], 'r', label='integrated eq 2')
    # plt.plot(log_time, [to_mag(lumi.Rayleign_Jeans(theta, (t * u.s).to(u.d))) for t in 10 ** log_time], 'b', label='blackbody')
    plt.legend()
    plt.grid()
    # plt.ylim(-12, -19)
    plt.gca().invert_yaxis()
    plt.show()


def test_figure_9_k():
    bands = ['U', 'B', 'V', 'R','I']
    colors = ['m', 'b', 'k', 'r', 'g']
    log_time = np.linspace(1.8, 6)
    for band, color in zip(bands, colors):
        print(band)
        lumi = FilteredL(band_filter=band)
        to_mag = lumi.luminosity_to_absolute_mag
        plt.plot(log_time, [to_mag(lumi.get_filtered_L(m12l15rot2, (t*u.s).to(u.d))) for t in 10**log_time], color, label=band)
        plt.ylim(-12.5, -18)
        plt.legend()
        plt.grid()
    plt.title('figure 9 Kozyreva')
    plt.xlabel('log time [log(s)]')
    plt.ylabel('magnitude [mags]')
    plt.show()

def test_fig_7():
    log_time = np.linspace(1.8, 6.2)
    lumi = FilteredL(band_filter='R', system='ABMag')
    to_mag = lumi.luminosity_to_absolute_mag
    eq_2 = lumi.eq_2
    nu = lumi.wave.length_to_frequency((length * 10 ** -9 * u.m).to_value(u.angstrom))
    print(nu)
    # plt.plot(log_time, [to_mag(lumi.get_filtered_L(theta, t*u.s)) for t in 10**log_time])
    # plt.plot(log_time, [to_mag(eq_2(theta, (t*u.s).to(u.d), nu)*nu) for t in 10**log_time])
    plt.plot(log_time, [to_mag(lumi.integrated_freq_L(theta, (t * u.s).to(u.d))) for t in 10 ** log_time])
    plt.ylim(-12, -20)
    plt.xlabel('log t [s]')
    plt.ylabel('M_AB')
    plt.title(f'Lambda = {length}nm')
    plt.grid()
    plt.show()


def test_random_stuff():
    a = 60 * u.s
    b = 1 * u.d
    print(b-a)

if __name__ == '__main__':
    test_fig_17_19_shusman(200)