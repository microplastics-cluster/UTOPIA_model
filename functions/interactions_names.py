import numpy as np
import pandas as pd

def extract_sinks(system_particle_object_list, SpeciesList):
    # Estimate Sinks (Losses and transformations and transport out of the box)

    loss_list_names = []
    loss_list_values= []
    
    for sp in system_particle_object_list:
        dict_sp=sp.RateConstants
        
        # replace none values
        for k in dict_sp.keys():
            if dict_sp[k] == None:
                dict_sp[k] = 0
            else:
                pass
        
        losses_names = []        
        for i in dict_sp.keys():
            if dict_sp[i]!= 0:
                losses_names.append(i)
            else:
                pass
        
        loss_list_names.append("+".join(losses_names))
        
        losses_values = []
        for i in dict_sp.keys():
            if i == "k_fragmentation":
                if type(dict_sp[i]) == tuple:
                    losses_values.append(dict_sp[i][0])
                else:
                    losses_values.append(dict_sp[i])
            elif i == "k_advective_transport" or i == "k_mixing":
                if type(dict_sp[i]) == tuple:
                    losses_values.append(sum(dict_sp[i]))
                else:
                    losses_values.append(dict_sp[i])
            else:
                losses_values.append(dict_sp[i])

        loss_list_values.append(-(sum(losses_values)))
    
    sinks=pd.DataFrame({"Loss_processes_names":loss_list_names,"loss_k_val_s-1":loss_list_values}, index=SpeciesList)
    
    
    return sinks


