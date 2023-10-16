"""Import modules and functions"""

import copy
import os
import pandas as pd

from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *
from functions.fillInteractions_df_fun_OOP import *
from functions.generate_modelObjects import*
from functions.generateRateConstants_particles_test import *
from functions.solver_SteadyState import*
from functions.extract_results import*
from functions.plot_results import*
from functions.massBalance import*

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

"""Generate objects"""

## choose input files to load
comp_impFile_name="\inputs_compartments_simple.csv"
comp_interactFile_name="\compartment_interactions.csv"
mp_imputFile_name="\inputs_microplastics.csv"

## generate objects
modelBoxes, system_particle_object_list, SpeciesList, process_inputs_df,spm, dict_comp,model_lists = generate_objects(inputs_path,boxName="Utopia", comp_impFile_name=comp_impFile_name, comp_interactFile_name=comp_interactFile_name, mp_imputFile_name=mp_imputFile_name)

"""Estimate rate constants per particle"""

## Estimate rate constants
for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)
    
## Create rate constants table:
RC_df= create_rateConstants_table(system_particle_object_list)

## Plot rate constants
plot_rate_constants(RC_df)

"""Build Matrix of interactions"""

interactions_df = fillInteractions_fun_OOP(system_particle_object_list, SpeciesList)

"""SOLVE SYSTEM OF ODES"""

# Choose input flow (in g per second or particle per second) 

q_mass_g_s=1
q_num_s=0


# Define particle imput (sp_imput):
    # Size fraction:
        # a= 0.5 um
        # b= 5 um
        # c= 50 um
        # d= 500 um
        # e= 5000 um
    # Aggregation state:
        # A= Free MP
        # B= heteroaggregatedMP
        # C= biofouled MP
        # D= biofouled and heteroaggregated MP
    # Compartment:
        # 0= 'Ocean_Surface_Water'
        # 1= 'Ocean_Mixed_Water'
        # 2= 'Ocean_Column_Water'
        # 3= 'Coast_Surface_Water'
        # 4= 'Coast_Column_Water'
        # 5= 'Surface_Freshwater'
        # 6= 'Bulk_Freshwater'
        # 7= 'Sediment'
        # 8= 'Urban_Soil_Surface'
        # 9= 'Urban_Soil'
        # 10= 'Background_Soil_Surface'
        # 11= 'Background_Soil'
        # 12= 'Agricultural_Soil_Surface'
        # 13= 'Agricultural_Soil'
        # 14= 'Air'
    
sp_imput="eA0_Utopia"

## Solve ODES

Results, R = solve_ODES_SS(system_particle_object_list=system_particle_object_list,q_mass_g_s=q_mass_g_s,q_num_s=q_num_s,sp_imput=sp_imput,interactions_df=interactions_df)

## Organize results

Results_comp_dict= extract_by_comp(R.reset_index(),model_lists["compartmentNames_list"])

Results_comp_organiced= extract_by_aggSt(Results_comp_dict,MPforms_list)

#Plot results in Total number of particle

particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}

for comp in Results_comp_organiced:
    plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])


## Solve mass balance and print result

massBalance(R,system_particle_object_list, q_mass_g_s)


###Modify rate constants by stablishing a time limit or chaging specific rate constant values using the change_RC_value function
# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

from functions.update_rateConstants import*

timeLimit_particles_RC(particle_object_list=system_particle_object_list,lim=0.000556,SpeciesList=SpeciesList,q_mass_g_s=q_mass_g_s,q_num_s=q_num_s,sp_imput=sp_imput,model_lists=model_lists,MPforms_list=MPforms_list)

change_RC_value(particle_object_list=system_particle_object_list,rc_name="k_sea_spreay_aerosol",rc_val=0,)