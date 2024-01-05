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


"""Define run parameters"""

## choose input files to load

comp_impFile_name = "\inputs_compartments.csv"
comp_interactFile_name = "\compartment_interactions.csv"
mp_imputFile_name = "\inputs_microplastics.csv"

boxName = "Utopia"

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

## Compartments:
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

# input flow (in g per second)
q_mass_g_s = 1

# particle imput
size_bin = "e"
compartment = "Air"
MP_form = "freeMP"
MP_density = "LowDensity"  # To be changed based on the MP imputs file

runName = (
    MP_density
    + "MP_Emissions_"
    + str(q_mass_g_s)
    + "g_s_"
    + MP_form
    + "_"
    + str(size_dict[size_bin])
    + "_nm_"
    + compartment
)

# Create directory with model run name under the current date directory where to save results

subDirectory = current_date + "_" + runName

path_run = os.path.join(path, subDirectory)

# check whether directory already exists
if not os.path.exists(path_run):
    os.mkdir(path_run)
    print("Folder %s created!" % path_run)
else:
    print("Folder %s already exists" % path_run)


"""Generate objects"""

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

surfComp_list = [c for c in dict_comp if "Surface" in c]

"""Estimate rate constants per particle"""

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)

#### Original results with no time limit

## create rate constants table:
RC_df = create_rateConstants_table(system_particle_object_list)
df4 = RC_df.fillna(0)

# Plot rate constants


"""(FIX RC for wet deposition, now its given as a list of rate constants per surface compartment only for dry deposition)This needs to be fixed also for the matrix of interactions and estimation of flows"""


# plot_rate_constants(df4)

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

# Plot and save RC heatmaps per compartment (Has to be fixed for the procesess wherer there are multiple values of rate constants)
for comp in dict_comp.keys():
    plot_heatmap_RC(comp, df4, path_RC)

"""Build Matrix of interactions"""

interactions_df = fillInteractions_fun_OOP(
    system_particle_object_list, SpeciesList, surfComp_list
)

from functions.fill_interactions_Knames import *

interactions_df_Knames = fillInteractions_Knames(
    system_particle_object_list, SpeciesList
)

"""SOLVE SYSTEM OF ODES"""

particle_compartmentCoding = dict(
    zip(
        model_lists["compartmentNames_list"],
        list(range(len(model_lists["compartmentNames_list"]))),
    )
)
comp_dict_inverse = {v: k for k, v in particle_compartmentCoding.items()}

sp_imput = (
    size_bin
    + particle_forms_coding[MP_form]
    + str(particle_compartmentCoding[compartment])
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
    + compartment
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

# Save results extended organised in excel sheets per compartment

results_extended_by_compartment_to_csv(path=path_run, results_dict=Results_comp_dict)

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
        if type(p.RateConstants[c]) == list:
            p.outFlow_mass_g_s[c] = [R * p.Pmass_g_SS for R in p.RateConstants[c]]
        else:
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
                for index, value in enumerate(surfComp_list):
                    if value == comp:
                        position = index
                df = tables_outputFlows[e_comp].loc[:, ["k_" + ele for ele in inpProc]]
                for proc in inpProc:
                    df["k_" + proc] = df["k_" + proc].apply(
                        lambda x: x[position] if isinstance(x, list) else x
                    )

                comp_input_flows.append(df)
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

comp_mass_balance = {}
for comp in list(dict_comp.keys()):
    comp_mass_balance[comp] = compartment_massBalance(
        comp=comp,
        tables_outputFlows=tables_outputFlows,
        PartMass_t0=PartMass_t0,
        comp_dict_inverse=comp_dict_inverse,
        dict_comp=dict_comp,
        tables_inputFlows=tables_inputFlows,
    )


""" Generate PDF report """  ## WORK IN PROGRESS
from functions.generate_pfd_report import *

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

create_pdf_report(df_list, figs, filename, text_elements)
