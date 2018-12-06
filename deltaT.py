###
#   Function to construct interpolator isotherm at a different temperature
#   Author: CoryMSimon@gmail.com
###
__author__ = 'Cory M. Simon'

import scipy.optimize
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IAST import InterpolatorIsotherm


def ConstructInterpolatorIsothermAtTnew(df, T0, Tnew, loading_key=None, pressure_key=None, hoa_key=None, fill_value=None):
    """
    Returns InterpolatorIsotherm from pure-component isotherm data measured at temperature T0 extrapolated to temperature Tnew using the heat of adsorption.

    :param df: DataFrame Pandas DataFrame with pure-component isotherm tabular data measured at temperature T0
    :param T0: Float temperature at which isotherm in df was measured (T0)
    :param Tf: Float temperature at which we desire to extrapolate the isotherm in df
    :param loading_key: String key for loading column in df
    :param pressure_key: String key for pressure column in df
    :param hoa_key: String key for heat of adsorption column in df
    :param fill_value: Float value of loading to assume when an attempt is made to interpolate at a pressure greater than the largest pressure observed in the data

    HOA needs to be in units kJ/mol.

    """
    if loading_key == None or pressure_key == None or hoa_key == None:
        raise Exception("Pass loading_key, hoa_key, and pressure_key, names of loading, heat of adsorption, and pressure cols in DataFrame.")
    if pressure_key == 'new_P':
        raise Exception("Change pressure column to something new plz")
    
    # for every point, shift pressures according to Classius-Clapyeron eqn
    R = 8.314 / 1000.0 # kJ/mol-K
    n = df.shape[0]
    df_new = pd.DataFrame()
    df_new[pressure_key] = np.zeros((n,))
    df_new[loading_key] = df[loading_key].values

    for i in range(n):
        df_new[pressure_key].iloc[i] = df[pressure_key].iloc[i] * np.exp(-df[hoa_key].iloc[i] / R * (1.0 / Tnew - 1.0 / T0))
    return InterpolatorIsotherm(df_new, loading_key=loading_key, pressure_key=pressure_key, fill_value=fill_value), df_new
