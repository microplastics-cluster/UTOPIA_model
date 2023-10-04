import copy
import os
import pandas as pd
from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *

from functions.generate_modelObjects import*

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

# Generate objects
## choose input files to load

modelBoxes, system_particle_object_list, SpeciesList, process_inputs_df,spm, dict_comp,model_lists = generate_objects(inputs_path,boxName="Utopia", comp_impFile_name="\inputs_compartments.csv", comp_interactFile_name="\compartment_interactions.csv", mp_imputFile_name="\inputs_microplastics.csv")


# Estimate rate constants per particle
from functions.generateRateConstants_particles import *

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)


###Modify rate constants by stablishing a time limit or chaging specific rate constant values using the change_RC_value function
# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

# particles_updated= change_RC_value(system_particle_object_list,rc_name,rc_val)

particles_updated = timeLimit_particles_RC(system_particle_object_list,0.000556)

RC_df_timeLim = create_rateConstants_table(particles_updated)

#Plot rate constants
plot_rate_constants(RC_df_timeLim)

# create rate constants table:

fileName = "rateConstantsUTOPIA_Test.csv"

rate_constants_df = create_rateConstants_table(particles_updated)

# Save rate contants dataframe as csv file

df4 = RC_df_timeLim.fillna(0)
df4.to_csv(fileName, index=False)


# Generate system of differentia equations (1-Matrix of interactions, 2-System of differential equations)

# Build Matrix of interactions
from functions.fillInteractions_df_fun_OOP import *

interactions_df = fillInteractions_fun_OOP(particles_updated, SpeciesList)


#Optional Check interactions dataframe by process:

from functions.fill_interactions_Knames import*


interactions_df_Knames=fillInteractions_Knames(
system_particle_object_list,SpeciesList
)


# """SOLVE SYSTEM OF ODES"""

from functions.solver_SteadyState import*
from functions.extract_results import*
from functions.plot_results import*

#Choose input flow (in g per second)
q_mass_g_s=1

Results, R = solve_ODES_SS(system_particle_object_list=particles_updated,q_mass_g_s=q_mass_g_s,q_num_s=0,sp_imput="eA0_Utopia",interactions_df=interactions_df)

from functions.massBalance import*
massBalance(R,system_particle_object_list, q_mass_g_s)

Results_comp_dict=extract_by_comp(R.reset_index(),model_lists["compartmentNames_list"])
Results_comp_organiced=extract_by_aggSt(Results_comp_dict,MPforms_list)

#Plot results
particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


for comp in Results_comp_organiced:
    plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])