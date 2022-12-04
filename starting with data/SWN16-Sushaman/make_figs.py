import emcee
from L_T_R import *
import pandas as pd
import matplotlib.pyplot as plt
import corner

def main():
    folder = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\SWN16-Sushaman\models\250 steps\no binning"
    model = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\SWN16-Sushaman\models\250 steps\no binning\after ext.h5"
    make_figs("", folder, model, 250)

def make_figs(data_path, folder, model, total_steps):
    print(model)
    reader = emcee.backends.HDFBackend(model, read_only=True)

    burnin = 0
    thin = 10
    flat_samples = reader.get_chain(discard=burnin, flat=True, thin=thin)


    labels = ["R500", "M15", "E51", "offset"]
    # for sheet_name in ['after ext binning', 'after ext']:
        # t, mag, mag_err = get_data(file_path=data_path, sheet_name=sheet_name)
        # get_lines_dots(flat_samples=flat_samples, x=t, y=mag, yerr=mag_err,
        #                file_path=folder + "\points_" + sheet_name)
        #
        # t = t[t < 1]
        # mag = mag[:len(t)]
        # mag_err = mag_err[:len(t)]
        #
        # get_lines_dots(flat_samples=flat_samples, n=n, x=t, y=mag, yerr=mag_err,
        #                file_path=folder + "\close_points_" + sheet_name, close_up=True)
    get_histograms(reader, file_path=folder+'\hist', labels=labels, total_steps=total_steps)
    get_corner(flat_samples, folder+'\corner', labels)


def get_data(file_path, sheet_name):
    data = pd.read_excel(file_path, sheet_name=sheet_name)

    t = np.array(data.loc[:, 'JD - 2457651.0[day]'])
    meas_mag = np.array(data.loc[:, 'V[mag]'])
    meas_mag_err = np.array(data.loc[:, 'error_V[mag]'])

    return t, meas_mag, meas_mag_err

def get_histograms(sampler, file_path, labels, total_steps):
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
        plt.title(f'binning {total_steps}')
    else:
        plt.title(f'no binning {total_steps}')
    plt.savefig(file_path)


# def get_lines_dots(flat_samples, x, y, yerr, n, file_path, start=0.01, end=3.5, close_up=False):
#     plt.clf()
#     inds = np.random.randint(len(flat_samples), size=100)
#     time = np.linspace(start, end)
#     for ind in inds:
#         theta = flat_samples[ind]
#         if close_up:
#             time = np.linspace(t_bigger_than(n, theta), 0.02)
#         v85, R13, M_e, offset = theta
#         plt.plot(time, get_mag(n, time + offset, theta), "C1", alpha=0.1)
#
#     plt.errorbar(x - offset, y, yerr=yerr, fmt=".k", capsize=0)
#
#     plt.xlabel(f'JD - {2457651.0+offset}[day]')
#     plt.ylabel('V[mag]')
#     if 'binning' in file_path:
#         plt.title(f'n={n} binning')
#     else:
#         plt.title(f'n={n} no binning')
#     plt.gca().invert_yaxis()
#     plt.grid()
#     plt.savefig(file_path)


def get_corner(flat_samples, file_path, labels):
    fig = corner.corner(flat_samples, labels=labels)
    fig.savefig(file_path)

if __name__ == '__main__':
    main()