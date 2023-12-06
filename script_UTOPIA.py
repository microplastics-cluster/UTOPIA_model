import copy
import os
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *
from functions.fillInteractions_df_fun_OOP import *
from functions.generate_modelObjects import *
from functions.generateRateConstants_particles import *
from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.plot_results import *
from functions.massBalance import *
from functions.fill_interactions_Knames import *


inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

outputs_path = os.path.join(os.path.dirname(__file__), "Results")

# get current date and time to store results
current_datetime = datetime.now().strftime("%Y-%m-%d %H-%M-%S")
current_date = datetime.now().strftime("%Y-%m-%d")

# Create directory with current date where to save results
directory = current_date
path = os.path.join(outputs_path, directory)

# check whether directory already exists
if not os.path.exists(path):
    os.mkdir(path)
    print("Folder %s created!" % path)
else:
    print("Folder %s already exists" % path)


"""Generate objects"""

## choose input files to load

comp_impFile_name = "\inputs_compartments.csv"
comp_interactFile_name = "\compartment_interactions.csv"
mp_imputFile_name = "\inputs_microplastics.csv"

boxName = "Utopia"

runName = "LowDensityMP"

# Create directory with model run name under the current date directory where to save results

subDirectory = current_date + "_" + runName

path_run = os.path.join(path, subDirectory)

# check whether directory already exists
if not os.path.exists(path_run):
    os.mkdir(path_run)
    print("Folder %s created!" % path_run)
else:
    print("Folder %s already exists" % path_run)


# Generate objects
(
    modelBoxes,
    system_particle_object_list,
    SpeciesList,
    process_inputs_df,
    spm,
    dict_comp,
    model_lists,
    particles_df,
    MPforms_list,
) = generate_objects(
    inputs_path,
    boxName=boxName,
    comp_impFile_name=comp_impFile_name,
    comp_interactFile_name=comp_interactFile_name,
    mp_imputFile_name=mp_imputFile_name,
)


"""Estimate rate constants per particle"""

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)

#### Original results with no time limit

## create rate constants table:
RC_df = create_rateConstants_table(system_particle_object_list)
df4 = RC_df.fillna(0)

# Plot rate constants
plot_rate_constants(df4)

# # Plot heatmaps of rate constants per compartment

# Create rate constants folder:
path_RC = os.path.join(path_run, "RateConstants")

# check whether directory already exists
if not os.path.exists(path_RC):
    os.mkdir(path_RC)
    print("Folder %s created!" % path_RC)
else:
    print("Folder %s already exists" % path_RC)

# Save rate constants dataframe

t_filename = os.path.join(path_RC, "RateConstants_table.csv")
df4.to_csv(t_filename, index=False)

# Plot and save RC heatmaps per compartment
for comp in dict_comp.keys():
    plot_heatmap_RC(comp, df4, path_RC)

"""Build Matrix of interactions"""

interactions_df = fillInteractions_fun_OOP(system_particle_object_list, SpeciesList)

from functions.fill_interactions_Knames import *

interactions_df_Knames = fillInteractions_Knames(
    system_particle_object_list, SpeciesList
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
size_dict = dict(zip(["a", "b", "c", "d", "e"], [0.5, 5, 50, 500, 5000]))
# Aggregation state:
# A= Free MP
# B= heteroaggregatedMP
# C= biofouled MP
# D= biofouled and heteroaggregated MP
particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))

MP_form_dict_reverse = {v: k for k, v in particle_forms_coding.items()}

# Compartment:
# 0= 'Ocean_Surface_Water'
# 1= 'Ocean_Mixed_Water'
# 2= 'Ocean_Column_Water'
# 3= 'Coast_Surface_Water'
# 4= 'Coast_Column_Water'
# 5= 'Surface_Freshwater'
# 6= 'Bulk_Freshwater'
# 7= 'Sediment_Freshwater'
# 8= 'Sediment_Ocean'
# 9= 'Sediment_Coast'
# 10= 'Urban_Soil_Surface'
# 11= 'Urban_Soil'
# 12= 'Background_Soil_Surface'
# 13= 'Background_Soil'
# 14= 'Agricultural_Soil_Surface'
# 15= 'Agricultural_Soil'
# 16= 'Air'
particle_compartmentCoding = dict(
    zip(
        model_lists["compartmentNames_list"],
        list(range(len(model_lists["compartmentNames_list"]))),
    )
)
comp_dict_inverse = {v: k for k, v in particle_compartmentCoding.items()}

# input flow (in g per second)
q_mass_g_s = 1

# particle imput
size_bin = "e"
comp = "Ocean_Surface_Water"
MP_form = "freeMP"

sp_imput = (
    size_bin
    + particle_forms_coding[MP_form]
    + str(particle_compartmentCoding[comp])
    + "_"
    + boxName
)

print(
    "Imput flow = "
    + str(q_mass_g_s)
    + " g/s of "
    + MP_form
    + " of size "
    + str(size_dict[size_bin])
    + " into the "
    + comp
    + " compartment"
)

R, PartMass_t0 = solve_ODES_SS(
    system_particle_object_list=system_particle_object_list,
    q_mass_g_s=q_mass_g_s,
    q_num_s=0,
    sp_imput=sp_imput,
    interactions_df=interactions_df,
)

# Reformat results (R) dataframe
R["Size_Fraction_um"] = [size_dict[x[0]] for x in R.index]
R["MP_Form"] = [MP_form_dict_reverse[x[1]] for x in R.index]
R["Compartment"] = [comp_dict_inverse[float(x[2:-7])] for x in R.index]

Results = R[
    [
        "Compartment",
        "MP_Form",
        "Size_Fraction_um",
        "mass_g",
        "number_of_particles",
        "concentration_g_m3",
        "concentration_num_m3",
    ]
]

# Solve mass balance and print result
massBalance(R, system_particle_object_list, q_mass_g_s)


# Test that there are no negative results
for i, idx in zip(R["mass_g"], R.index):
    if i < 0:
        print("negative values in the solution for " + idx)
    else:
        pass

# Estimate mass and number fractions and extract ranking tables of the species with higest fractions to understand the distribution of the particles in the system by mass and number of particles

Results_extended, mf_shorted, nf_shorted = estimate_fractions(Results)


# Organise results in dictionary for plotting

Results_comp_dict = extract_by_comp(
    Results_extended.reset_index(), particle_compartmentCoding
)
Results_comp_organiced = extract_by_aggSt(Results_comp_dict, particle_forms_coding)

# Plot results in Total number of particles and Total mass

# Create Steady State results folders:
path_SteadyState_mass = os.path.join(path_run, "SteadyState_mass_distribution")
path_SteadyState_number = os.path.join(path_run, "SteadyState_number_distribution")


# check whether directory already exists
if not os.path.exists(path_SteadyState_mass):
    os.mkdir(path_SteadyState_mass)
    print("Folder %s created!" % path_SteadyState_mass)
else:
    print("Folder %s already exists" % path_SteadyState_mass)

if not os.path.exists(path_SteadyState_number):
    os.mkdir(path_SteadyState_number)
    print("Folder %s created!" % path_SteadyState_number)
else:
    print("Folder %s already exists" % path_SteadyState_number)


for comp in Results_comp_organiced:
    plot_bySize_total_number_particles(
        Results_comp_organiced,
        comp,
        model_lists["dict_size_coding"],
        path=path_SteadyState_number,
    )
    plot_bySize_total_mass(
        results_dict=Results_comp_organiced,
        comp_name=comp,
        dict_size_coding=model_lists["dict_size_coding"],
        path=path_SteadyState_mass,
    )
    # plot_by(
    #     results_dict=Results_comp_organiced,
    #     comp_name=comp,
    #     dict_size_coding=model_lists["dict_size_coding"],
    #     plot_by="concentration_g_m3",
    # )


print("Distribution of mass in the system")
print(mf_shorted[:10])
df_massDistribution = mf_shorted[:10]

# Save the table of mass distribution
massDitribution_filename = os.path.join(
    path_SteadyState_mass, "SS_mass_distribution.csv"
)
df_massDistribution.to_csv(massDitribution_filename)


print("distribution of particle number in the system")
print(nf_shorted[:10])
df_numberDistribution = nf_shorted[:10]

# Save the table of number distribution
numberDitribution_filename = os.path.join(
    path_SteadyState_number, "SS_number_distribution.csv"
)
df_numberDistribution.to_csv(numberDitribution_filename)


# Mass distribution by compartment
mass_frac_100 = []
num_frac_100 = []
for comp in list(dict_comp.keys()):
    mass_frac_100.append(
        sum(Results_extended[Results_extended["Compartment"] == comp]["mass_fraction"])
        * 100
    )
    num_frac_100.append(
        sum(
            Results_extended[Results_extended["Compartment"] == comp]["number_fraction"]
        )
        * 100
    )

mass_dist_comp = pd.DataFrame(columns=["Compartments"])
mass_dist_comp["Compartments"] = list(dict_comp.keys())
mass_dist_comp["%_mass"] = mass_frac_100
mass_dist_comp["%_number"] = num_frac_100


# Plot Mass distribution by compartment

sns.barplot(data=mass_dist_comp, x="Compartments", y="%_mass").set(
    title="Mass distribution"
)
plt.ylabel("% of total mass")
plt.yscale("log")
plt.ylim(0.01, 100)
plt.xticks(rotation=90)
# Save the plot
massDistribPlot_filename = os.path.join(
    path_SteadyState_mass, "mass_by_compartment.png"
)
plt.savefig(massDistribPlot_filename, bbox_inches="tight")
plt.show()
plt.close()

# Plot Number distribution by compartment
sns.barplot(data=mass_dist_comp, x="Compartments", y="%_number").set(
    title="Particle number distribution"
)
plt.ylabel("% of total number")
plt.yscale("log")
plt.ylim(0.01, 100)
plt.xticks(rotation=90)
# Save the plot
numberDistribPlot_filename = os.path.join(
    path_SteadyState_number, "number_by_compartment.png"
)
plt.savefig(numberDistribPlot_filename, bbox_inches="tight")
plt.show()
plt.close()


### MASS BALANCE PER COMPARTMENT###

# Estimate mass flows due to the different particle fate process (transfer between compartments, elimination and transformation processes)

# Estimate outflows
for p in system_particle_object_list:
    p.outFlow_mass_g_s = {}
    for c in p.RateConstants:
        p.outFlow_mass_g_s[c] = p.RateConstants[c] * p.Pmass_g_SS

# Tables of output flows per compartmet
tables_outputFlows = {}
for c in list(dict_comp.keys()):
    part_dic = {}
    for p in system_particle_object_list:
        if p.Pcompartment.Cname == c:
            part_dic[p.Pcode] = pd.DataFrame.from_dict(
                p.outFlow_mass_g_s, orient="index"
            )
    tables_outputFlows[c] = pd.concat(part_dic, axis=1).transpose()
for k in tables_outputFlows:
    tables_outputFlows[k] = (
        tables_outputFlows[k].reset_index(level=1).drop("level_1", axis=1)
    )

# Estimate imput flows from transport from other compartments

##Tables of recieving flows through transport from other compartments
tables_inputFlows = {}
for comp in list(dict_comp.keys()):
    comp_input_flows = []
    for e_comp in dict_comp:
        if comp in dict_comp[e_comp].connexions:
            inpProc = dict_comp[e_comp].connexions[comp]
            if type(inpProc) == list:
                comp_input_flows.append(
                    tables_outputFlows[e_comp].loc[:, ["k_" + ele for ele in inpProc]]
                )
            else:
                comp_input_flows.append(
                    tables_outputFlows[e_comp].loc[:, "k_" + inpProc].to_frame()
                )
        else:
            pass

    tables_inputFlows[comp] = pd.concat(comp_input_flows).fillna(0)


# Create folder for saving ouput and input mass flows
path_mass_flows = os.path.join(path_run, "Compartment_mass_flows")

# check whether directory already exists
if not os.path.exists(path_mass_flows):
    os.mkdir(path_mass_flows)
    print("Folder %s created!" % path_mass_flows)
else:
    print("Folder %s already exists" % path_mass_flows)

for comp in tables_outputFlows:
    flows_tables_to_csv(
        comp,
        tables_outputFlows,
        tables_inputFlows,
        MP_form_dict_reverse,
        size_dict,
        path=path_mass_flows,
    )


## Compartment mass balance --> Check mass balance function (when done manually it works.... adding inputs from table of inputs + emissions - outflows (selecting the rigth processess))


for comp in list(dict_comp.keys()):
    compartment_massBalance(
        comp=comp,
        tables_outputFlows=tables_outputFlows,
        PartMass_t0=PartMass_t0,
        comp_dict_inverse=comp_dict_inverse,
        dict_comp=dict_comp,
        tables_inputFlows=tables_inputFlows,
    )


""" Generate PDF report """

filename = runName + "_" + current_date
text_elements = {
    "plastic_density_kg_m3": system_particle_object_list[0].Pdensity_kg_m3,
    "imput_flow_g_s": q_mass_g_s,
    "particle_emissions_form": MP_form,
    "particle_emissions_size_nm": size_dict[size_bin],
    "recieving_compartment": comp,
}
df_list = [df_massDistribution, df_numberDistribution, df4]
figs = ["rateConstants.png", "massDistribution.png", "numberDistribution.png"]


""" Save results """

# Create directory with current date
directory = current_date

path = os.path.join(outputs_path, directory)

# check whether directory already exists
if not os.path.exists(path):
    os.mkdir(path)
    print("Folder %s created!" % path)
else:
    print("Folder %s already exists" % path)

# Save results in directory

## Particle properties
particles_df.to_csv(
    os.path.join(path, "particle_properties_" + runName + "_" + current_date + ".csv")
)

## Process inputs
process_inputs_df.to_csv(
    os.path.join(path, "process_inputs_" + runName + "_" + current_date + ".csv")
)

## Rate constants
df4.to_csv(
    os.path.join(path, "rate_constants_" + runName + "_" + current_date + ".csv")
)

## Interactions matrix
# interactions_df.to_csv(
#     os.path.join(path, "interactions_matrix_" + runName + "_" + current_date + ".csv")
# )
# interactions_df_Knames.to_csv(
#     os.path.join(
#         path, "interactions_matrix_Knames_" + runName + "_" + current_date + ".csv"
#     )
# )

## Results at steady state

Results_extended.to_csv(
    os.path.join(path, "Results_SS_" + runName + "_" + current_date + ".csv")
)

# Save figures?

# for comp in mass_conc_figures.keys():
#     mass_conc_figures[comp].saveFigure()


###Modify rate constants by stablishing a time limit or chaging specific rate constant values using the change_RC_value function
# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

# particles_updated= change_RC_value(system_particle_object_list,#rc_name="k_sediment_resuspension",rc_val=1E-7)

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

### Code to check

# Estimate mass flows of particles out of a compartment due to transfer to other compartments (rate constants of transfer * C at steady steat)
for p in system_particle_object_list:
    p.outFlow_mass_g_s = {}
    for c in p.Pcompartment.connexions.values():
        if type(c) != list:
            p.outFlow_mass_g_s[c] = p.RateConstants["k_" + c] * p.Pmass_g_SS
        else:
            for e in c:
                p.outFlow_mass_g_s[e] = p.RateConstants["k_" + e] * p.Pmass_g_SS

# Print tables of output flows per compartmet
tables_outputFlows = {}
for c in list(dict_comp.keys()):
    comp_out_flow_table = pd.DataFrame(
        columns=["MP_form", "MP_size"] + list(dict_comp[c].connexions.values())
    )
    for p in system_particle_object_list:
        if p.Pcompartment.Cname == c:
            comp_out_flow_table.loc[len(comp_out_flow_table)] = [
                p.Pform,
                p.Pcode[0],
            ] + [p.outFlow_mass_g_s[z] for z in list(dict_comp[c].connexions.values())]
    tables_outputFlows[c] = comp_out_flow_table
