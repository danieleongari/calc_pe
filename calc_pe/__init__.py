from __future__ import absolute_import
import sys
import logging
import pyiast
import numpy as np
import pandas as pd
from six.moves import range

# pylint: disable=too-many-arguments,too-many-locals,too-many-statements


def wcvf(y, T, P, vf, ms):
    """ Returns the moles of CO2 per kg of adsorbent in the void fraction (wcvf)
    considering the ideal gas law

    :param y: CO2 fraction in the gas phase (-)
    :param T: temperature of the system (K)
    :param P: total pressure (Pa)
    :param vf: void fraction of the bed
    :param ms: mass of adsorbent in per volume of bed (kg/m3)
    """
    R = 8.3144621  # J/mol/K
    return P * y * vf / (R * T * ms)


def ConstructInterpolatorIsothermAtTnew(df,
                                        T_0,
                                        T_new,
                                        loading_key=None,
                                        pressure_key=None,
                                        hoa_key=None,
                                        fill_value=None):
    """
    Returns InterpolatorIsotherm from pure-component isotherm data measured at
    temperature T0 extrapolated to temperature Tnew using the heat of adsorption.

    :param df: DataFrame Pandas DataFrame with pure-component isotherm tabular
               data measured at temperature T0
    :param T_0: Float temperature at which isotherm in df was measured (T0)
    :param T_new: Float temperature at which we desire to extrapolate the isotherm
               in df
    :param loading_key: String key for loading column in df
    :param pressure_key: String key for pressure column in df
    :param hoa_key: String key for heat of adsorption column in df
    :param fill_value: Float value of loading to assume when an attempt is made
                       to interpolate at a pressure greater than the largest
                       pressure observed in the data

    HOA needs to be in units kJ/mol.
    """

    R = 8.3144621 / 1000.0  # kJ/mol/K
    n = df.shape[0]  # number of values in the isotherm
    df_new = pd.DataFrame()
    df_new[pressure_key] = np.zeros((n, ))
    df_new[loading_key] = df[loading_key].values

    # for every point, shift pressures according to Classius-Clapyeron eqn
    for i in range(n):
        p_0 = df[pressure_key].iloc[i]  # Pa
        HoA = df[hoa_key].iloc[i]  # kJ/mol (negative)
        p_new = p_0 * np.exp(HoA / R * (1.0 / T_new - 1.0 / T_0))
        df_new[pressure_key].iloc[i] = p_new

    # return a pyiast.isotherm object
    return pyiast.InterpolatorIsotherm(df_new,
                                       loading_key=loading_key,
                                       pressure_key=pressure_key,
                                       fill_value=fill_value), df_new


def totalQ(Ta, Td, Pa, Pd, ya, yd, iso_Ta, iso_Td, cp, vf, ms, qa, wcvf_a,
           iso_df):
    """ Returns [0] Q_{thermal} computed from eq. 1 in 10.1039/c4ee02636e (kJ/mol)
              [1] Qs, sensible heat (kJ/mol)
              [2] Qd, enthalpy of desorption (kJ/mol, defined as positive)
              [3] WCv, mass CO2 per m3 of bed (kg/m3)
              [4] pur, final CO2 purity (-)

  """
    CO2_KGMOL = 0.044  # kg/mol
    # extract or compute the adsorption conditions if not already stored
    if not qa:
        logging.debug(
            'Ta = {:.3e}, Pa = {:.3e}, ya = {:.3e}, yd = {:.3e}'.format(
                Ta, Pa, ya, yd))
        wcvf_a['CO_2'] = wcvf(ya, Ta, Pa, vf, ms)
        wcvf_a['N_2'] = wcvf(1 - ya, Ta, Pa, vf, ms)
        # gas uptake @ adsorption
        qa['CO_2'], qa['N_2'] = pyiast.iast(
            np.array([ya, 1 - ya]) * Pa,
            [iso_Ta['CO_2'], iso_Ta['N_2']],
            verboseflag=False,
            warningoff=True,  # no complain about extrapolation
        )
        logging.debug("qa['CO_2']: {:.3e}, qa['N_2']: {:.3e}".format(
            qa['CO_2'], qa['N_2']))
    # Compute desorption conditions
    wcvf_d = {}
    wcvf_d['CO_2'] = wcvf(yd, Td, Pd, vf, ms)
    wcvf_d['N_2'] = wcvf(1 - yd, Td, Pd, vf, ms)
    qd = {}
    # gas uptake @ desorption
    qd['CO_2'], qd['N_2'] = pyiast.iast(
        np.array([yd, 1 - yd]) * Pd,
        [iso_Td['CO_2'], iso_Td['N_2']],
        verboseflag=False,
        warningoff=True,  # no complain about extrapolation
        adsorbed_mole_fraction_guess=[0.99999999, 0.00000001])
    logging.debug("qd['CO_2']: {:.3e}, qd['N_2']: {:.3e}".format(
        qd['CO_2'], qd['N_2']))
    wc, wct = {}, {}
    for m in ['CO_2', 'N_2']:
        wc[m] = qa[m] - qd[
            m]  # working capacity in the adsorbent (mol_component/kg_ADS)
        wct[m] = wc[m] + wcvf_a[m] - wcvf_d[
            m]  # working capacity total = adsorbent + void
    WCv = wct['CO_2'] * ms * CO2_KGMOL  # kgCO2/m3
    logging.debug("wct['CO_2'] = {:.3e} => WCv = {:.3e}".format(
        wct['CO_2'], WCv))

    if WCv < 0:
        logging.debug("NEGATIVE wct['CO_2']")
        return np.nan, np.nan, np.nan, WCv, np.nan

    # Compute volumetric working capacity, and final purity
    pur = wct['CO_2'] / (wct['CO_2'] + wct['N_2'])
    logging.debug('WCv = {:.3e}, pur = {:.3e}'.format(WCv, pur))
    # Compute Q_{thermal} (Qt) = sensible_heat (Qs) + desorption_enthalpy (Qd)
    Qs = cp * ms * (Td - Ta) / WCv  # J/kgCO2
    Qd = 0
    for m in ['CO_2', 'N_2']:
        # Get the enthalpy from the input data, at the closest conditions to desorption (assuming constant HoA with loading)
        idx_closest_Pd = (iso_df[m]['pressure(Pa)'] - Pd).abs().argsort()[0]
        h = -iso_df[m].loc[idx_closest_Pd][
            'HoA(kJ/mol)']  # positive number, we have to provide this enthalpy
        Qd += h * 1e3 * wc[m] * ms / WCv  # J/kgCO2
    Qt = Qs + Qd  # J/kgCO2
    logging.debug('Qs = {:.3e}, Qd = {:.3e}, Qt = {:.3e}'.format(Qs, Qd, Qt))
    return Qt, Qs, Qd, WCv, pur


def totalWcomp(Pd, pur):
    """ Compression work (J/kgCO2) from P_desorption (Pa) to 150bar. This uses
    a 'functional representation' as described in 10.1039/c4ee02636e.
    """
    plog = -5.4433e5 / np.log(1e6 / 1e3)
    Wp = 7.1723e5 + np.log(Pd / 1000) * plog
    mp = (2.313 - 2.102) / (np.sqrt(1e6) - np.sqrt(1e3))
    mpr = 1.102 + mp * (np.sqrt(Pd) - np.sqrt(1e3))
    Wcomp = Wp * (1 + mpr * (1. / pur - 1))
    return Wcomp


def totalPE(Pa, Pd, Ta, Td, ya, yd, iso_Ta, iso_Td, cp, ms, vf, eleff, qa,
            wcvf_a, iso_df):
    """ Total PE, computed from eq. 2 in 10.1039/c4ee02636e (eleff = 'carnot') """
    Qt, Qs, Qd, WCv, pur = totalQ(Ta, Td, Pa, Pd, ya, yd, iso_Ta, iso_Td, cp,
                                  vf, ms, qa, wcvf_a, iso_df)
    if WCv < 0:  # Negative working capacity: return NaN for all values but WCv
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, WCv, np.nan
    # Compute compression work from Pd to 150bar
    Wcomp = totalWcomp(Pd, pur)  # J/kgCO2
    logging.debug('Wcomp = {:.3e}'.format(Wcomp))
    if eleff == 'carnot':
        nu = (Td + 10 - 283) / (Td + 10)
        logging.debug('nu = {:.3e}'.format(nu))
        Qteff = 0.75 * Qt * nu
        Qseff = 0.75 * Qs * nu
        Qdeff = 0.75 * Qd * nu
    elif eleff == 'linear':
        Qteff = Qt * (-0.08037 + 0.002326 * (Td - 273 + 10))
        Qseff = Qs * (-0.08037 + 0.002326 * (Td - 273 + 10))
        Qdeff = Qd * (-0.08037 + 0.002326 * (Td - 273 + 10))
    else:
        sys.exit('WARNING: eleff parameter is not in [carnot,linear]')
    logging.debug('Qt = {:.3e}, Qteff = {:.3e}'.format(Qt, Qteff))
    PE = Qteff + Wcomp  # J/kgCO2
    return PE, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur


def mainPE(gasin, rho, vf, process, cp, yd, eleff, opt, T_iso, iso_df):  # noqa
    """Main script to compute CO2 parasitic energy from single component
    CO2 and N2 isotherms

    :str gasin: choice of the input mixture
    :float rho: density of the adsorbent (kg/m3)
    :float vf:  void fraction of the bed, e.g., 0.35
    :str process: choice of the process, e.g., TPSA
    :float cp: heat capacity of the adsorbent (J/kg/K), e.g., 985.0
    :float yd: target CO2 purity, e.g., 0.99
    :str eleff: choice of the method to compute electr. eff., e.g., 'carnot'
    :str opt: choice of the parameter to optimize for the process, e.g., 'PE'
    :{'CO_2': float,'N_2': float} T_iso: temperatures of the used isotherms
    :{'CO_2': df,'N_2': df} iso_df: dataframes with isotherms, see examples
    """

    # Initialize calculation
    ms = rho * (1. - vf)  # Mass (kg) of adsorbent in per m3 of bed

    if gasin == 'coal':
        totE = 6631.2 / 1.80  # kJ/kg_CO2 (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.14  # 14:86 molar ratio
        Ta = 313.0
        Pa = 101325.0
    elif gasin == 'ng':
        totE = 21023.26 / 3.23  # kJ/kg_CO2 (10.1039/C4EE02636E, SI, pag. 6)
        ya = 0.04  # 4% molar CO2
        Ta = 313.0
        Pa = 101325.0
    elif gasin == 'air':
        totE = np.nan  # ill defined
        ya = 0.0004  # 400ppm molar CO2
        Ta = 288.0
        Pa = 101325.0
    iso_Ta = {}
    for m in ['CO_2', 'N_2']:
        # Extrapolate the pure-component isotherm at Ta
        iso_Ta[m], _ = ConstructInterpolatorIsothermAtTnew(
            iso_df[m],
            T_iso[m],
            Ta,
            loading_key='loading(mol/kg)',
            pressure_key='pressure(Pa)',
            hoa_key='HoA(kJ/mol)',
            fill_value=iso_df[m]['loading(mol/kg)'].max())

    # Set the range and the step for the Td and Pd, tested to find the min PE:
    # the range was changed from {10.1039/c4ee02636e, pag. 4136, 2nd col},
    # to have less points and be faster.
    Td_min = Ta + 20.0
    Td_step = 10.0
    Td_max = Td_min + 100.0
    Pd_min = 0.01 * 101325.0
    Pd_step = 0.02 * 101325.0
    Pd_max = 1.5 * 101325.0
    if process == 'PSA':
        Td_range = np.array([Td_min])
        Pd_range = np.arange(Pd_min, Pd_max, Pd_step)
    elif process == 'TSA':
        Td_range = np.arange(Td_min, Td_max, Td_step)
        Pd_range = np.array([Pa])
    elif process == 'TPSA':
        Td_range = np.arange(Td_min, Td_max, Td_step)
        Pd_range = np.arange(Pd_min, Pd_max, Pd_step)
    else:
        sys.exit('WARNING: process parameter is not in [PSA,TSA,TPSA]')
    # Compute the PE at different Td, Pd
    data = []  # collect all the data at different Pd, Td
    qa = {}
    wcvf_a = {}
    iso_Td = {}
    for Td in Td_range:
        for m in ['CO_2', 'N_2']:
            # Extrapolate the pure-component isotherm at Td
            iso_Td[m], _ = ConstructInterpolatorIsothermAtTnew(
                iso_df[m],
                T_iso[m],
                Td,
                loading_key='loading(mol/kg)',
                pressure_key='pressure(Pa)',
                hoa_key='HoA(kJ/mol)',
                fill_value=iso_df[m]['loading(mol/kg)'].max())
        for Pd in Pd_range:
            # Compute the PE @ Td and Pd
            logging.debug('**** Evaluating: Td={:.1f}, Pd={:.3f} *****'.format(
                Td, Pd))
            PE, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur = totalPE(
                Pa, Pd, Ta, Td, ya, yd, iso_Ta, iso_Td, cp, ms, vf, eleff, qa,
                wcvf_a, iso_df)
            if WCv > 0:  #Positive working capacity
                logging.debug('PE = {:.3e}'.format(PE))
                data.append(
                    [PE, Pd, Td, Qteff, Qt, Qseff, Qdeff, Wcomp, WCv, pur])

    # Find the optimal conditions required and combine the results in a dict
    if not data:  # All wc are negative
        results_dict = {'process_feasible': False}
        return results_dict
    data = np.array(data)
    if opt == 'PE':
        data_opt = data[np.argmin(data.T[0])]
    elif opt == 'Q':
        data_opt = data[np.argmin(data.T[3])]
    elif opt == 'WC':
        data_opt = data[np.argmax(data.T[8])]
    elif opt == 'pur':
        data_opt = data[np.argmax(data.T[9])]
    else:
        sys.exit('WARNING: opt parameter is not in [PE,Q,WC,pur]')
    logging.debug('data_opt:')
    logging.debug(data_opt)
    results_dict = {
        'process_feasible': True,
        'PE': data_opt[0] / 1e6,
        'PE_units': 'MJ/kg',
        'PE_descr': 'parassitic energy at optimal process conditions',
        'Pd': data_opt[1] / 101325.0,
        'Pd_units': 'bar',
        'Pd_descr': 'desorption pressure',
        'Td': data_opt[2],
        'Td_units': 'K',
        'Td_descr': 'desorption temperature',
        'eloss': data_opt[0] / 1e3 / totE,
        'eloss_units': 'kJ/kJ',
        'eloss_descr': 'fraction of electricity loss (np.nan for air)',
        'Qteff': data_opt[3] / 1e6,
        'Qteff_units': 'MJ/kgCO2',
        'Qteff_descr': 'heat requirement',
        'Qt': data_opt[6] / 1e6,
        'Qt_units': 'MJ/kgCO2',
        'Qt_descr': 'Not-effective heat requirement',
        'Qs': data_opt[4] / 1e6,
        'Qs_units': 'MJ/kgCO2',
        'Qs_descr': 'Sensible heat',
        'Qd': data_opt[5] / 1e6,
        'Qd_units': 'MJ/kgCO2',
        'Qd_descr': 'desorption enthalpy',
        'Wcomp': data_opt[7] / 1e6,
        'Wcomp_units': 'MJ/kgCO2',
        'Wcomp_descr': 'compression work',
        'WCv': data_opt[8],
        'WCv_units': 'kgCO2/m3',
        'WCv_descr': 'volumetric working capacity',
        'WCg': data_opt[8] / ms,
        'WCg_units': 'kgCO2/kg',
        'WCg_descr': 'gravimetric working capacity',
        'Pur': data_opt[9],
        'Pur_units': 'mol/mol',
        'Pur_descr': 'fraction of CO2 purity',
    }
    return results_dict
