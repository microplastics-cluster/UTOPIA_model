from functions.fillInteractions_df_fun_OOP import eliminationProcesses
from helpers.helpers import num_to_mass

import pandas as pd
import numpy as np

    
def massBalance(R,system_particle_object_list, q_mass_g_s):
    
    # Estimate looses: loss processess=[discorporation, burial]
    elimination_rates = []
    for p in system_particle_object_list:
        elimination_rates.append(p.RateConstants["k_discorporation"])
        
    #mass at Steady state
    m_ss=R["mass_g"]   
    
    #output flow
    out_low_g_s=sum(elimination_rates*m_ss)
    
    print("Difference inflow-outflow = "+ str(q_mass_g_s-out_low_g_s))
    
        
  