from functions import *

def test():
    data_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\excel files\combined_data.xlsx"

    burnin = 200
    thin = 20

    labels = ["v85", "R13", "M_e", "offset"]
    for sheet_name in ['after ext binning', 'after ext']:
        t, mag, mag_err = get_data(file_path=data_path, sheet_name=sheet_name)

        t = t[t < 1]
        print(t)
        mag = mag[:len(t)]
        mag_err = mag_err[:len(t)]




def make_figs(data_path, model, folder, n):
    reader = emcee.backends.HDFBackend(folder+'/'+model)

    burnin = 600
    thin = 60
    flat_samples = reader.get_chain(discard=burnin, flat=True, thin=thin)

    labels = ["v85", "R13", "M_e", "offset"]
    for sheet_name in ['after ext binning', 'after ext']:
        t, mag, mag_err = get_data(file_path=data_path, sheet_name=sheet_name)
        get_lines_dots(flat_samples=flat_samples, n=n, x=t, y=mag, yerr=mag_err,
                       file_path=folder + "\points_" + sheet_name)

        t = t[t < 1]
        mag = mag[:len(t)]
        mag_err = mag_err[:len(t)]

        get_lines_dots(flat_samples=flat_samples, n=n, x=t, y=mag, yerr=mag_err,
                       file_path=folder + "\close_points_" + sheet_name, close_up=True)

    get_histograms(reader, file_path=folder+'\hist_'+sheet_name, labels=labels, n=n)
    get_corner(flat_samples, folder+'\corner_'+sheet_name, labels, n)

def main():
    data_path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\excel files\combined_data.xlsx"
    new_models = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\SW16\models\with t bigger"
    for folder in os.listdir(new_models):
        folder = new_models + '/' + folder
        model = os.listdir(folder)[0]
        # print(model)
        if '3.0' in model:
            n = 3
        else:
            n = 1.5
        make_figs(data_path, model, folder, n)

if __name__ == '__main__':
    main()

