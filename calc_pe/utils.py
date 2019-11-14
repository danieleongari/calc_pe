from __future__ import absolute_import
from __future__ import print_function
import logging


def get_PE_results_string(adsorbent_name, res):
    """Combine PE results into a single line string."""
    if res['process_feasible']:
        results_str = '{:s}: '.format(adsorbent_name)
        results_str += 'PE(MJ/kg)= {:.3f} '.format(res['PE'])
        results_str += 'Pd(bar)= {:.2f} '.format(res['Pd'])
        results_str += 'Td(K)= {:.1f} '.format(res['Td'])
        results_str += 'EL(J/J)= {:.3f} '.format(res['eloss'])
        results_str += 'Q(MJ/kg)= {:.3f} '.format(res['Qteff'])
        results_str += 'Wcomp(MJ/kg)= {:.3f} '.format(res['Wcomp'])
        results_str += 'WCv(kg/m3)= {:.3f} '.format(res['WCv'])
        results_str += 'WCg(kg/kg)= {:.3f} '.format(res['WCg'])
        results_str += 'pur(mol/mol)= {:.3f}'.format(res['Pur'])
    else:
        results_str = '{:s}: Unfeasible process!'.format(adsorbent_name)
    # Print in the log file
    logging.debug(results_str)

    return results_str


def printPE(adsorbent_name, res):
    """Print PE_results_string on screen."""
    results_str = get_PE_results_string(adsorbent_name, res)
    print(results_str)
