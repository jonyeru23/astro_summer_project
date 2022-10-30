import pytest
from L_T_R import *


#############################################################################
@pytest.fixture
def wave_factory():
    return {
        'V': Wave(),
        'U': Wave('U'),
        'R': Wave('R')
    }


def test_length_to_freq(wave_factory):
    assert abs((wave_factory['V'].central_freq / 1e10) - 5.4601e+14 / 1e10) < 1
    assert abs((wave_factory['U'].central_freq / 1e10) - 8.3309e+14 / 1e10) < 1
    assert abs((wave_factory['R'].central_freq / 1e10) - 4.6123e+14 / 1e10) < 1


def test_range_freq(wave_factory):
    for wave in wave_factory.values():
        # print(type(wave.ob.wave[0]))
        freqs = wave.get_range_freq(10)
        for i in range(9):
            assert freqs[i + 1] - freqs[i] > 0

        assert wave.ob.wave == np.flip(wave.frequency_to_length(freqs))


###########################################################################################################

# @pytest.fixture
# def luminosity_factory():
#     return {
#         'V': FilteredL(100),
#         'U': FilteredL(100, 'U', 'ABMag'),
#         'R': FilteredL(100, 'R'),
#         'not_filtered': L()
#     }

@pytest.fixture
def progenitors_factory():
    m12l15rot2 = [812 / 500, 8.09 / 15, 1.1, 0]
    m12l5rot2 = [345 / 500, 9.15 / 15, 1.1, 0]
    m15l15rot0 = [1024 / 500, 11.33 / 15, 1.1, 0]
    m15l5rot8 = [268 / 500, 6.1 / 15, 1.1, 0]
    return {
        'm12l15rot2': m12l15rot2,
        'm12l5rot2': m12l5rot2,
        'm15l15rot0': m15l15rot0,
        'm15l5rot8': m15l5rot8
    }


# def test_mag_to_lum(luminosity_factory):
#     for band in ['V', 'U']:
#         mag_to_flux = luminosity_factory[band].mag_to_flux
#         luminosity_to_mag = luminosity_factory[band].luminosity_to_absolute_mag
#         assert abs(mag_to_flux(luminosity_to_mag(1)) - 1) < 0.0001
#
#
# def test_right_mag(luminosity_factory, progenitor_factory):
#     pass

test_range_freq(wave_factory)
