from get_model import *

def main():
    # sheet_name = input('Enter sheet name:')
    # steps = int(input('Enter steps:'))
    steps = 100
    path = r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\excel files\combined_data.xlsx"
    for sheet_name in ['after ext binning', 'after ext']:
        sampler = Sampler(path, sheet_name=sheet_name)
        sampler.write_sampler(f"{sheet_name}.h5", steps=steps)


if __name__ == '__main__':
    main()

