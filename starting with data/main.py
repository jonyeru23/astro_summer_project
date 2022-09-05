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
    t, mag, mag_err = get_data(file_path='excel files/combined_data.xlsx', sheet_name='after ext binning')
    for n in [3/2, 3]:
        get_sampler(f"n={n}.h5", n=n, nwalkers=16, steps=4000, x=t, y=mag, yerr=mag_err)




if __name__ == '__main__':
    main()