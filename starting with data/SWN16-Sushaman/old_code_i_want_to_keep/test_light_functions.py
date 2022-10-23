import pandas as pd
import pytest
from light_functions import *

@pytest.fixture
def web_data(**kwargs):
    def _data_factory(file_name):
        return pd.read_csv(file_name)
    return _data_factory


def test_L_obs(web_data):
    theta = [R_x(1, 500), M_x(1, 15), E_x(1, 51), 0]
    data = web_data(r"C:\Users\User\OneDrive - mail.tau.ac.il\Desktop\אוניברסיטה\אסטרו נודר\פרויקט קיץ\התחלה של קוד\astro_summer_project\starting with data\SWN16-Sushaman\test_data\R=1_M=1_E=1.csv")
    for test_l, t in zip(np.array(data.loc[:, 'L_bol[erg/s]']), np.array(data.loc[:, 't[s]'])):
        print(test_l - L_obs(theta, (t*u.s).to(u.d)))
        assert abs(test_l - L_obs(theta, (t*u.s).to(u.d))) < 10**3
