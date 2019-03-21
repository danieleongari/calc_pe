import logging

def printPE(adsorbent_name, res):
    """ Combine PE results into a single line string and print them on screen """
    if res['process_feasible']:
        results_str="{:s}: ".format(adsorbent_name)
        results_str+="PE(MJ/kg)= {:.3f} ".format(res['PE'])
        results_str+="Pd(bar)= {:.2f} ".format(res['Pd'])
        results_str+="Td(K)= {:.1f} ".format(res['Td'])
        results_str+="EL(J/J)= {:.3f} ".format(res['eloss'])
        results_str+="Q(MJ/kg)= {:.3f} ".format(res['Qteff'])
        results_str+="Wcomp(MJ/kg)= {:.3f} ".format(res['Wcomp'])
        results_str+="WCv(kg/m3)= {:.3f} ".format(res['WCv'])
        results_str+="WCg(kg/kg)= {:.3f} ".format(res['WCg'])
        results_str+="pur(mol/mol)= {:.3f}".format(res['Pur'])
    else:
        results_str="{:s}: Unfeasible process!".format(adsorbent_name)
    # Print in the log file
    logging.debug(results_str)
    # Print on screen
    print(results_str)
