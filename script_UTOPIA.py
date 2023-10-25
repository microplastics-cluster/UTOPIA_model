import copy
import os
import pandas as pd
from datetime import datetime

from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *
from functions.fillInteractions_df_fun_OOP import *
from functions.generate_modelObjects import*
from functions.generateRateConstants_particles import *
from functions.solver_SteadyState import*
from functions.extract_results import*
from functions.plot_results import*
from functions.massBalance import*

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

outputs_path = os.path.join(os.path.dirname(__file__), "Results")

# get current date and time to store results
current_datetime = datetime.now().strftime("%Y-%m-%d %H-%M-%S")
current_date = datetime.now().strftime("%Y-%m-%d")


"""Generate objects"""

## choose input files to load

comp_impFile_name = "\inputs_compartments.csv"
comp_interactFile_name="\compartment_interactions.csv"
mp_imputFile_name="\inputs_microplastics.csv"

boxName="Utopia"

runName="LowDensityMP"

modelBoxes, system_particle_object_list, SpeciesList, process_inputs_df,spm, dict_comp,model_lists,particles_df,MPforms_list = generate_objects(inputs_path,boxName=boxName,comp_impFile_name=comp_impFile_name, comp_interactFile_name=comp_interactFile_name, mp_imputFile_name=mp_imputFile_name)


"""Estimate rate constants per particle"""

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)

#### Original results with no time limit

## create rate constants table:
RC_df= create_rateConstants_table(system_particle_object_list)
df4 = RC_df.fillna(0)

# Plot rate constants
plot_rate_constants(df4)

"""Build Matrix of interactions"""

interactions_df = fillInteractions_fun_OOP(system_particle_object_list, SpeciesList)

from functions.fill_interactions_Knames import*

interactions_df_Knames=fillInteractions_Knames(system_particle_object_list,SpeciesList
)

"""SOLVE SYSTEM OF ODES"""

# Choose input flow (in g per second) 
# Define particle imput (sp_imput):

# Size fraction:
    # a= 0.5 um
    # b= 5 um
    # c= 50 um
    # d= 500 um
    # e= 5000 um
size_dict=dict(zip(["a","b","c","d","e"],[0.5,5,50,500,5000]))
# Aggregation state:
    # A= Free MP
    # B= heteroaggregatedMP
    # C= biofouled MP
    # D= biofouled and heteroaggregated MP
particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))

MP_form_dict_reverse={v: k for k, v in particle_forms_coding.items()}

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
particle_compartmentCoding = dict(
    zip(model_lists["compartmentNames_list"], list(range(len(model_lists["compartmentNames_list"]))))
)
comp_dict_inverse = {v: k for k, v in particle_compartmentCoding.items()}

# input flow (in g per second)
q_mass_g_s=1

# particle imput
size_bin= "e"
comp="Ocean_Surface_Water"
MP_form="freeMP"

sp_imput=size_bin+particle_forms_coding[MP_form]+str(particle_compartmentCoding[comp])+"_"+boxName

print("Imput flow = "+str(q_mass_g_s) +" g/s of "+ MP_form + " into the "+ comp +" compartment" )

R = solve_ODES_SS(system_particle_object_list=system_particle_object_list,q_mass_g_s=q_mass_g_s,q_num_s=0,sp_imput=sp_imput,interactions_df=interactions_df)

#Reformat results (R) dataframe
R["Size_Fraction_um"]=[size_dict[x[0]] for x in R.index]
R["MP_Form"]=[MP_form_dict_reverse[x[1]] for x in R.index]
R["Compartment"]=[comp_dict_inverse[float(x[2:-7])] for x in R.index]

Results=R[['Compartment','MP_Form','Size_Fraction_um','mass_g', 'number_of_particles', 'concentration_g_m3',
       'concentration_num_m3']]

#Solve mass balance and print result
massBalance(R,system_particle_object_list, q_mass_g_s)


Results_comp_dict=extract_by_comp(R.reset_index(),particle_compartmentCoding)
Results_comp_organiced=extract_by_aggSt(Results_comp_dict,particle_forms_coding)

# Test that there are no negative results
for i,idx in zip(R["mass_g"],R.index):
    if i<0:
        print("negative values in the solution for " + idx)
    else:
        pass


#Plot results in Total number of particles and Total mass
for comp in Results_comp_organiced:
    plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])
    plot_bySize_total_mass(results_dict=Results_comp_organiced,comp_name=comp,dict_size_coding=model_lists["dict_size_coding"])


# Estimate mass and number fractions and extract ranking table of the species with higest fractions to understan the distribution of the particles in the system

Results_extended=estimate_fractions(Results)

    
    
""" Save results """

# Create directory with current date
directory = current_date

path= os.path.join(outputs_path, directory) 

# check whether directory already exists
if not os.path.exists(path):
  os.mkdir(path)
  print("Folder %s created!" % path)
else:
  print("Folder %s already exists" % path)
  
# Save results in directory

## Particle properties
particles_df.to_csv(os.path.join(path,"particle_properties_"+runName+"_"+current_date+".csv"))

## Process inputs
process_inputs_df.to_csv(os.path.join(path,"process_inputs_"+runName+"_"+current_date+".csv"))

## Rate constants
df4.to_csv(os.path.join(path,"rate_constants_"+runName+"_"+current_date+".csv"))

## Interactions matrix
interactions_df.to_csv(os.path.join(path,"interactions_matrix_"+runName+"_"+current_date+".csv"))
interactions_df_Knames.to_csv(os.path.join(path,"interactions_matrix_Knames_"+runName+"_"+current_date+".csv"))

## Results at steady state

Results.to_csv(os.path.join(path,"Results_SS_"+runName+"_"+current_date+".csv"))
    



###Modify rate constants by stablishing a time limit or chaging specific rate constant values using the change_RC_value function
# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

#particles_updated= change_RC_value(system_particle_object_list,#rc_name="k_sediment_resuspension",rc_val=1E-7)

# particles_updated = timeLimit_particles_RC(system_particle_object_list,0.000556)


# RC_df_timeLim = create_rateConstants_table(particles_updated)

# #Plot rate constants
# plot_rate_constants(RC_df_timeLim)

# # create rate constants table:

# fileName = "rateConstantsUTOPIA_Test.csv"

# # Save rate contants dataframe as csv file

# df4 = RC_df_timeLim.fillna(0)
# #df4.to_csv(fileName, index=False)


# # Generate system of differentia equations (1-Matrix of interactions, 2-System of differential equations)

# # Build Matrix of interactions

# interactions_df = fillInteractions_fun_OOP(particles_updated, SpeciesList)


# #Optional Check interactions dataframe by process:

# from functions.fill_interactions_Knames import*


# #interactions_df_Knames=fillInteractions_Knames(
# #system_particle_object_list,SpeciesList
# #)


# # """SOLVE SYSTEM OF ODES"""

# #Choose input flow (in g per second)
# q_mass_g_s=1
# sp_imput="eA0_Utopia"

# R = solve_ODES_SS(system_particle_object_list=particles_updated,q_mass_g_s=q_mass_g_s,q_num_s=0,sp_imput=sp_imput,interactions_df=interactions_df)

# from functions.massBalance import*
# massBalance(R,system_particle_object_list, q_mass_g_s)

# Results_comp_dict=extract_by_comp(R.reset_index(),model_lists["compartmentNames_list"])
# Results_comp_organiced=extract_by_aggSt(Results_comp_dict,MPforms_list)

# #Plot results
# particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


# for comp in Results_comp_organiced:
#     plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])