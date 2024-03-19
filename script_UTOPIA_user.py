import copy
import os
from datetime import datetime

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *
from functions.fillInteractions_df_fun import *
from functions.generate_modelObjects import *
from functions.generateRateConstants_particles import *
from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.plot_results import *
from functions.massBalance import *
from functions.fill_interactions_Knames import *
from functions.exposure_indicators_calculation import *
from functions.generate_MPinputs_table import *
from functions.save_results import *


inputs_path = os.path.join(os.path.dirname(__file__), "inputs")


"""Define run parameters"""

## Define microplastics physical properties

# The user can also select a preloaded file instead of typing in the values. In this case the user wont need to run the code between lines 29 and 34 and neither the code between lines 42 and 50. The user will have to run line 56 with the selected input file

MPdensity_kg_m3 = 1580
MP_composition = "PVC"
shape = "sphere"  # Fixed for now
N_sizeBins = 5  # Fixed, should not be changed. The 5 size bins are generated as being one order of magnitude appart and cover the range from mm to nm(i.e. 5000um, 500um, 50um, 5um, 0.5um)
big_bin_diameter_um = 5000  # This size can not be bigger than 10 mm (10000um) or smaller than 1 mm(1000um)
runName = MP_composition

# write microplastics inputs file
mp_imputFile_name = write_MPinputs_table(
    MPdensity_kg_m3,
    MP_composition,
    shape,
    N_sizeBins,
    big_bin_diameter_um,
    runName,
    inputs_path,
)

## Environmental Characteristics

## Suspended particulates properties

spm_diameter_um = 0.5
spm_density_kg_m3 = 2000


## choose input files to load

comp_impFile_name = "\inputs_compartments.csv"  # Preloaded values, the user should be able to create its own inputs_compartments.csv file (via donwloading the file and typing news values without chaing the structure of the file) when a new file wants to be used the name should be changed here
comp_interactFile_name = (
    "\compartment_interactions.csv"  # Fixed, should not be modified
)
# mp_imputFile_name = os.path.join(inputs_path, "inputs_microplastics.csv") #Choose one existing input file to load

boxName = "Utopia"  # fixed, do not modify

"""Generate objects"""

# Generate objects
(
    system_particle_object_list,
    SpeciesList,
    spm,
    dict_comp,
    model_lists,
    particles_df,
) = generate_objects(
    inputs_path,
    boxName=boxName,
    MPforms_list=MPforms_list,
    comp_impFile_name=comp_impFile_name,
    comp_interactFile_name=comp_interactFile_name,
    mp_imputFile_name=mp_imputFile_name,
    spm_diameter_um=spm_diameter_um,
    spm_density_kg_m3=spm_density_kg_m3,
)

surfComp_list = [c for c in dict_comp if "Surface" in c]

## Microplastics weathering properties

## Select fragmentation style
"""estimate fragmentation relation between size bins using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html). Each particle fractions into fragments of smaller sizes and the distribution is expresses via the fragment size distribution matrix fsd. # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class """

if N_sizeBins == 5:
    frag_styles_dict = {
        "sequential_fragmentation": np.array(
            [
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
            ]
        ),
        "erosive_fragmentation": np.array(
            [
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [0.99, 0.01, 0, 0, 0],
                [0.999, 0, 0.001, 0, 0],
                [0.9999, 0, 0, 0.0001, 0],
            ]
        ),
        "mixed_fragmentation": np.array(
            [
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [0.5, 0.5, 0, 0, 0],
                [0.6, 0.2, 0.2, 0, 0],
                [0.7, 0.15, 0.1, 0.05, 0],
            ]
        ),
        "no_fragmentation": np.array(
            [
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ]
        ),
    }

    frag_style = "mixed_fragmentation"

    fsd = frag_styles_dict[frag_style]
    sizes = [list(model_lists["dict_size_coding"].keys())]
    fsd_df = pd.DataFrame(fsd, index=sizes, columns=sizes)

    # Save the fsd matrix
    fsd_filename = os.path.join(inputs_path, "fsd.csv")
    fsd_df.to_csv(fsd_filename)


else:
    print(
        "Fragmetation size distribution not defined for this number of size fractions, please define manually the fsd matrix via the fsd.csv file"
    )
    fsd_df = pd.read_csv(os.path.join(inputs_path, "fsd.csv"), index_col=0)
    fsd = fsd_df.to_numpy()

# optionally the user can type its own fsd matrix following the desciption above


## Weahering processes input parameters

# Generate the process inputs table based on the given model structure (created model boxes, compartments and particles)

## Degradation half time: thalf_deg_d
"Values used in Domercq et al. 2021, go to publication for more details on the selection of these values and asumptions made"
# Assumptions:
# Heteroaggregated particles degrade 10 times slower than the free MPs
# Biofouled particles degrade 5 times slower than the free MPs

# thalf_deg_d_dict = {
#     "freeMP": 5000,
#     "heterMP": 50000,
#     "biofMP": 25000,
#     "heterBiofMP": 100000,
# } #default values

# # Save the fsd matrix
# t_half_deg_filename = os.path.join(inputs_path, "t_half_deg.csv")
# t_half_deg_df = pd.DataFrame(list(thalf_deg_d_dict.items()), columns=['MP_form', 'thalf_deg_d'])
# t_half_deg_df.to_csv(t_half_deg_filename,index=False)


# If user wants to modify the default thalf_deg_d_dict, they can do so here or through the csv file and upload it

# Read the CSV file into a DataFrame
t_half_deg_filename = os.path.join(inputs_path, "t_half_deg.csv")
t_half_deg_df = pd.read_csv(t_half_deg_filename)

# Convert the DataFrame to a dictionary
thalf_deg_d_dict = t_half_deg_df.set_index("MP_form")["thalf_deg_d"].to_dict()

# Heteroaggregation attachment efficiency: alpha_heter.
alpha_heter_filename = os.path.join(inputs_path, "alpha_heter.csv")
alpha_heter_df = pd.read_csv(alpha_heter_filename)
alpha_hetr_dict = alpha_heter_df.set_index("MP_form")["alpha_heter"].to_dict()

# Timescale for fragmentation of the biggest size fraction (mp5): tfrag_gen_d

# t_frag_gen_df=

process_inputs_df = create_inputsTable_UTOPIA(
    inputs_path, model_lists, thalf_deg_d_dict, alpha_hetr_dict
)

"""Revisit create inputs table function...assumptions to be discussed and parameters to be added"""

## Emission Scenario

# Choose input flow (in g per second)
# Define particle imput (sp_imput): the user has to define in wich form and size the particles are released into the environment and specify the input flow for each compartment

# Size fraction:
# for the preloaded scenario:
# a= 0.5 um
# b= 5 um
# c= 50 um
# d= 500 um
# e= 5000 um
import string

size_codes = [letter for letter in string.ascii_lowercase[0:N_sizeBins]]
size_dict = dict(zip(size_codes, model_lists["dict_size_coding"].values()))

size_bin = "e"  # Chosse from size_dict


# Aggregation state (MP form):
# A= Free MP
# B= heteroaggregatedMP
# C= biofouled MP
# D= biofouled and heteroaggregated MP
MPforms_list = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))
MP_form_dict_reverse = {v: k for k, v in particle_forms_coding.items()}

MP_form = "freeMP"  # Choose from MPforms_list above

# input flow (in g per second) for each compartment the User should specify here the input flows per compartment
q_mass_g_s_dict = {
    "Ocean_Surface_Water": 1,
    "Ocean_Mixed_Water": 0,
    "Ocean_Column_Water": 0,
    "Coast_Surface_Water": 0,
    "Coast_Column_Water": 0,
    "Surface_Freshwater": 0,
    "Bulk_Freshwater": 0,
    "Sediment_Freshwater": 0,
    "Sediment_Ocean": 0,
    "Sediment_Coast": 0,
    "Urban_Soil_Surface": 0,
    "Urban_Soil": 0,
    "Background_Soil_Surface": 0,
    "Background_Soil": 0,
    "Agricultural_Soil_Surface": 0,
    "Agricultural_Soil": 0,
    "Air": 0,
}

input_flow_filename = os.path.join(inputs_path, "inputFlows.csv")
input_flows_df = pd.DataFrame(
    list(q_mass_g_s_dict.items()), columns=["compartment", "q_mass_g_s"]
)
input_flows_df.to_csv(input_flow_filename, index=False)

# input_flows_df = pd.read_csv(input_flow_filename)
# q_mass_g_s_dict=input_flows_df.set_index('compartment')['q_mass_g_s'].to_dict()

saveName = (
    MP_composition
    + "_MP_Emissions_"
    + MP_form
    + "_"
    + str(size_dict[size_bin])
    + "_nm_"
)

"""Estimate rate constants per particle"""

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp, fsd)


## create rate constants table:
RC_df = create_rateConstants_table(system_particle_object_list)
df4 = RC_df.fillna(0)

# Plot rate constants

"""(FIX RC for wet deposition, now its given as a list of rate constants per surface compartment only for dry deposition and wet depossition is turned off)This needs to be fixed also for the matrix of interactions and estimation of flows"""


"""Build Matrix of interactions"""

interactions_df = fillInteractions_fun_OOP(
    system_particle_object_list, SpeciesList, surfComp_list
)


"""SOLVE SYSTEM OF ODES"""

particle_compartmentCoding = dict(
    zip(
        model_lists["compartmentNames_list"],
        list(range(len(model_lists["compartmentNames_list"]))),
    )
)
comp_dict_inverse = {v: k for k, v in particle_compartmentCoding.items()}

sp_imputs = []
q_mass_g_s = []
for compartment in q_mass_g_s_dict.keys():

    sp_imputs.append(
        size_bin
        + particle_forms_coding[MP_form]
        + str(particle_compartmentCoding[compartment])
        + "_"
        + boxName
    )
    q_mass_g_s.append(q_mass_g_s_dict[compartment])

imput_flows_g_s = dict(zip(sp_imputs, q_mass_g_s))


R, PartMass_t0 = solve_ODES_SS(
    system_particle_object_list=system_particle_object_list,
    q_num_s=0,
    imput_flows_g_s=imput_flows_g_s,
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

# Total number of particles and Total mass

print("Distribution of mass in the system")
print(mf_shorted[:10])
df_massDistribution = mf_shorted[:10]


print("distribution of particle number in the system")
print(nf_shorted[:10])
df_numberDistribution = nf_shorted[:10]

# Mass distribution by compartment
mass_frac_100 = []
num_frac_100 = []
mass_conc_g_m3 = []
num_conc = []
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
    mass_conc_g_m3.append(
        sum(
            Results_extended[Results_extended["Compartment"] == comp][
                "concentration_g_m3"
            ]
        )
    )
    num_conc.append(
        sum(
            Results_extended[Results_extended["Compartment"] == comp][
                "concentration_num_m3"
            ]
        )
    )

mass_dist_comp = pd.DataFrame(columns=["Compartments"])
mass_dist_comp["Compartments"] = list(dict_comp.keys())
mass_dist_comp["%_mass"] = mass_frac_100
mass_dist_comp["%_number"] = num_frac_100
mass_dist_comp["Concentration_g_m3"] = mass_conc_g_m3
mass_dist_comp["Concentration_num_m3"] = num_conc


### MASS BALANCE PER COMPARTMENT###

# Estimate mass flows due to the different particle fate process (transfer between compartments, elimination and transformation processes)

from functions.generate_compartmentFlows_tables import *

# Estimate outflows
tables_outputFlows = estimate_outFlows(system_particle_object_list, dict_comp)


# Estimate imput flows from transport from other compartments
tables_inputFlows = estimate_inFlows(tables_outputFlows, dict_comp, surfComp_list)


## Compartment mass balance

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

# Print compartment mass balance table
comp_mass_balance_df = pd.DataFrame.from_dict(comp_mass_balance, orient="index")
print(comp_mass_balance_df)

comp_mass_balance_df["Mass balance"] = [
    comp_mass_balance_df["Inflow"][c] - comp_mass_balance_df["Outflow"][c]
    for c in comp_mass_balance_df.index
]

# Add total steady state mass and number of particles concentrations to dataframe

# comp_mass_balance_df["Total Mass (g)"] = [sum(Results_comp_dict[c].mass_g) for c in comp_mass_balance_df.index]
# comp_mass_balance_df["Total Number of Particles"] = [sum(Results_comp_dict[c].number_of_particles) for c in comp_mass_balance_df.index]
comp_mass_balance_df["Concentration (g/m3)"] = [
    sum(Results_comp_dict[c].concentration_g_m3) for c in comp_mass_balance_df.index
]
comp_mass_balance_df["Concentration (N/m3)"] = [
    sum(Results_comp_dict[c].concentration_num_m3) for c in comp_mass_balance_df.index
]


""" Estimate exposure indicators """

# Overall residence time
overall_residence_time_calculation(tables_outputFlows, Results_extended)


# Save results

outputs_path = os.path.join(os.path.dirname(__file__), "Results")

# Create directory with current date where to save results

# get current date and time to store results
current_date = datetime.now().strftime("%Y-%m-%d")
directory = current_date
path = os.path.join(outputs_path, directory)

# Create directory with model run name under the current date directory where to save results

subDirectory = current_date + "_" + saveName

path_run = os.path.join(path, subDirectory)

store_results(
    path,
    outputs_path,
    saveName,
    path_run,
    df4,
    Results_comp_dict,
    Results_comp_organiced,
    model_lists,
    df_massDistribution,
    df_numberDistribution,
    mass_dist_comp,
    tables_outputFlows,
    tables_inputFlows,
    MP_form_dict_reverse,
    size_dict,
    comp_mass_balance_df,
)

""" Generate PDF report """  ## WORK IN PROGRESS
# from functions.generate_pfd_report import *

# filename = saveName + "_" + current_date
# text_elements = {
#     "plastic_density_kg_m3": system_particle_object_list[0].Pdensity_kg_m3,
#     "imput_flow_g_s": q_mass_g_s,
#     "particle_emissions_form": MP_form,
#     "particle_emissions_size_nm": size_dict[size_bin],
#     "recieving_compartment": comp,
# }
# df_list = [df_massDistribution, df_numberDistribution, df4]
# figs = ["rateConstants.png", "massDistribution.png", "numberDistribution.png"]

# create_pdf_report(df_list, figs, filename, text_elements)
