from functions import *

def test():
    # labels = ["v85", "R13", "M_e", "offset"]
    t, mag, mag_err = get_data(file_path='excel files/combined_data.xlsx', sheet_name='after ext binning')
    # sampler = get_sampler(3/2, 16, 10, t, mag, mag_err)
    # get_histograms(sampler, "hist.png", labels=labels)
    # flat_samples = sampler.get_chain(flat=True)
    # draw_lines_dots(flat_samples, x=t, y=mag, yerr=mag_err, n=3/2, file_name='test_lines.png')
    # get_corner(flat_samples=flat_samples, file_name='test_corner.png', labels=labels)
    # sampler = get_sampler('test.h5', 3/2, 16, 10, t, mag, mag_err)

def main():
    # for n in [3/2, 3]:
    #     for sheet_name in ['after ext binning', 'after ext']:
    sheet_name = input('Enter sheet name:')
    n = int(input('Enter n:'))
    steps = int(input('Enter steps:'))
    t, mag, mag_err = get_data(file_path='combined_data.xlsx', sheet_name=sheet_name)
    get_sampler(f"n={n}_{sheet_name}.h5", n=n, nwalkers=16, steps=steps, x=t, y=mag, yerr=mag_err)






if __name__ == '__main__':
    main()