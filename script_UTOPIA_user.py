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
from functions.exposure_indicators_calculation import *
from functions.generate_MPinputs_table import *
from functions.save_results import *
from functions.loop_CTD_calculation import *
from functions.generate_compartmentFlows_tables import *
from functions.emission_fractions_calculation import *
from helpers.helpers import *
from functions.fillInteractions_dictionaries import *

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")


"""Define run parameters"""

## Define microplastics physical properties

# The user can also select a preloaded file instead of typing in the values. In this case the user wont need to run the code between lines 29 and 34 and neither the code between lines 42 and 50. The user will have to run line 56 with the selected input file

MPdensity_kg_m3 = 980
MP_composition = "PE"
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

# From Kooi et al. (2017)
v_a = 2.0e-16  # Volume of 1 algal cell [m-3]
r_a = ((3.0 / 4.0) * (v_a / math.pi)) ** (1.0 / 3.0)  # radius of algae [m]

spm_radius_um = r_a * 1e6
spm_density_kg_m3 = 1388  # REF: Kooi et al. (2017)


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
    spm_radius_um=spm_radius_um,
    spm_density_kg_m3=spm_density_kg_m3,
)

surfComp_list = [c for c in dict_comp if "Surface" in c]

## Microplastics weathering properties

## Select fragmentation style
"""estimate fragmentation relation between size bins using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html). Each particle fractions into fragments of smaller sizes and the distribution is expresses via the fragment size distribution matrix fsd. # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class """

# We provide a slider (to be done) where the user can select a fragmentation style by means of choosing a value of FI (fragmentation index) between 0 and 1 that describes two scenarios :

#     - Erosive fragmentation (FI=0): In this scenario the particles are being eroded on their surface and therefore most of their mass remain in their same size fraction and samall fraction in going to the samllest size bins. Its representative fsd is:

#         [[0, 0, 0, 0, 0],

#         [1, 0, 0, 0, 0],

#         [0.99, 0.01, 0, 0, 0],

#         [0.999, 0, 0.001, 0, 0],

#         [0.9999, 0, 0, 0.0001, 0],]

#     - Sequential fragmentation (FI=1): in this scenario each size fraction breacks down completely into the next smallest size bin.
#     Its representative fsd is:

#         [[0, 0, 0, 0, 0],

#         [1, 0, 0, 0, 0],

#         [0, 1, 0, 0, 0],

#         [0, 0, 1, 0, 0],

#         [0, 0, 0, 1, 0],]

#     By choosing a value between 0 and 1 the user can select a fragmentation style in between both extremes. (i.e. FI=0.5 will represent the mixed fragmentation style)

## Select a value for FI in the range 0-1 from a slider where we shoud indicate Erosive Fragmentation under the value 0 and Sequential Fragmentation under the value 1, like in the following dictionary:
# frag_styles_dict = {0:"erosive_fragmentation",0.5:"mixed_fragmentation",1:"sequential_fragmentation"}

FI = 0.5  # from slider or user imput

# Generate the fsd matrix
fsd = generate_fsd_matrix(FI)
# Create a dataframe from the fsd matrix
sizes = [list(model_lists["dict_size_coding"].keys())]
fsd_df = pd.DataFrame(fsd, index=sizes, columns=sizes)

# Save the fsd matrix (not sure if we need to save it-- to be revisited)
fsd_filename = os.path.join(inputs_path, "fsd.csv")
fsd_df.to_csv(fsd_filename)


## Weathering processes input parameters

# Generate the process inputs table based on the given model structure (created model boxes, compartments and particles)

##### Define Weathering processes input parameters

##### Degradation half time: thalf_deg_d

# The assumptions made for the definition of these degradation times: (NEW assumptions)
#     - Heteroaggregated particles degrade 10 times slower than the free MPs
#     - Biofouled particles degrade 2 times faster than the free MPs
#     - Both degradation and fragmentation rates are compartment dependent : we assume that in the surface water compartments both degradation and fragmentation are fastest, in soil surface and deeper water compartments both rates are 10 times slower (factor_deepWater_soilSurface) and in sediments and deeper soil compartments they both are 100 times slower (factor_sediment)

t_half_deg_free = 66000  # in days (10 times slower than the rate of degradation (to form dissolved organics) shown in Pfohl et al. 2023 for TPU-arom)
heter_deg_factor = 10
biof_deg_factor = 1 / 2

t_half_deg_heter = t_half_deg_free * heter_deg_factor
t_half_deg_biof = t_half_deg_free * biof_deg_factor
t_half_deg_biofHeter = t_half_deg_free * biof_deg_factor * heter_deg_factor

thalf_deg_d_dict = {
    "freeMP": t_half_deg_free,
    "heterMP": t_half_deg_heter,
    "biofMP": t_half_deg_biof,
    "heterBiofMP": t_half_deg_biofHeter,
}

factor_deepWater_soilSurface = 10
factor_sediment = 100

# t_half_deg_filename = os.path.join(inputs_path, "t_half_deg.csv")
# t_half_deg_df = pd.DataFrame(list(thalf_deg_d_dict.items()), columns=['MP_form', 'thalf_deg_d'])
# t_half_deg_df.to_csv(t_half_deg_filename,index=False)

# Heteroaggregation attachment efficiency: alpha_heter.
alpha_heter_filename = os.path.join(inputs_path, "alpha_heter.csv")
alpha_heter_df = pd.read_csv(alpha_heter_filename)
alpha_hetr_dict = alpha_heter_df.set_index("MP_form")["alpha_heter"].to_dict()

# Timescale for fragmentation

# The fragmentation timescales are deteremined from the stablished fragmentation half time of 36.5 days for the biggest size fraction in free form in the surface water compartments following the parameters chosen in Domercq et al. 2021.

# In UTOPIA we include fragmentation of the heteroaggregated MPs as being 100 slower than fragmentation of the Free MPs and breackup of biofouled and heteroaggregated will be two times slowed of those only heteroaggregated, following the same assumption as for free and biofouled. These values are used in the Domercq et al. 2021 paper and they are asumptions made from lack of current knowlegde

t_frag_gen_FreeSurfaceWater = 36.5  # in days
biof_frag_factor = 2
heter_frag_factor = 100


process_inputs_df = create_inputsTable_UTOPIA(
    inputs_path,
    model_lists,
    thalf_deg_d_dict,
    alpha_hetr_dict,
    t_frag_gen_FreeSurfaceWater,
    biof_frag_factor,
    heter_frag_factor,
    factor_deepWater_soilSurface,
    factor_sediment,
    save_op="save",
)

"""Revisit create inputs table function...assumptions to be discussed and parameters to be added"""

## Emission Scenario

# Choose input flow (in g per second)
# Define particle imput (sp_imput): the user has to define in wich form and size the particles are released into the environment and specify the input flow for each compartment

# Size fraction:
# for the preloaded scenario:
# a= 0.5 um =mp1
# b= 5 um
# c= 50 um
# d= 500 um
# e= 5000 um =mp5
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

input_flow_g_s = 1

emiss_comp = "Ocean_Surface_Water"

q_mass_g_s_dict = {
    "Ocean_Surface_Water": 0,
    "Ocean_Mixed_Water": 0,
    "Ocean_Column_Water": 0,
    "Coast_Surface_Water": 0,
    "Coast_Column_Water": 0,
    "Surface_Freshwater": 0,
    "Bulk_Freshwater": 0,
    "Sediment_Freshwater": 0,
    "Sediment_Ocean": 0,
    "Sediment_Coast": 0,
    "Beaches_Soil_Surface": 0,
    "Beaches_Deep_Soil": 0,
    "Background_Soil_Surface": 0,
    "Background_Soil": 0,
    "Impacted_Soil_Surface": 0,
    "Impacted_Soil": 0,
    "Air": 0,
}

q_mass_g_s_dict[emiss_comp] = input_flow_g_s

# If emissions are made to several compartments type the input flows corresponding to the specific compartments in the following dictionary:

# q_mass_g_s_dict = {
#     "Ocean_Surface_Water": 0,
#     "Ocean_Mixed_Water": 0,
#     "Ocean_Column_Water": 0,
#     "Coast_Surface_Water": 0,
#     "Coast_Column_Water": 0,
#     "Surface_Freshwater": 0,
#     "Bulk_Freshwater": 0,
#     "Sediment_Freshwater": 0,
#     "Sediment_Ocean": 0,
#     "Sediment_Coast": 0,
#     "Urban_Soil_Surface": 0,
#     "Urban_Soil": 0,
#     "Background_Soil_Surface": 0,
#     "Background_Soil": 0,
#     "Agricultural_Soil_Surface": 0,
#     "Agricultural_Soil": 0,
#     "Air": 0,
# }

# If inputs are made to different size classess and MP forms a new dictionary has to be used (TO be done)

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
    + "_FI_"
    + str(FI)
)

# Print model run summary

print("Model run: ")
print("Emissions flow (g/s): ", input_flow_g_s)
desired_key = next(key for key, value in q_mass_g_s_dict.items() if value > 0)
print("Receiving compartment/s: ", desired_key)
print("Emitted MP density (kg/m3): ", MPdensity_kg_m3)
print("Emitted MP shape: ", shape)
print("Emitted MP form: ", MP_form)
print("Emitted MP size (um): ", size_dict[size_bin])
print(saveName)


"""Estimate rate constants per particle"""

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp, fsd, process_inputs_df)


## create rate constants table:
RC_df = create_rateConstants_table(system_particle_object_list)
df4 = RC_df.fillna(0)

# Plot rate constants (not implemented anymore)

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

q_num_s = [
    mass_to_num(v, p.Pvolume_m3, p.Pdensity_kg_m3) if v != 0 else 0
    for k, v in zip(imput_flows_g_s.keys(), imput_flows_g_s.values())
    for p in system_particle_object_list
    if k == p.Pcode
]

# imput_flows_num_s = dict(zip(sp_imputs, q_num_s))


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


# Mass distribution by compartment
results_by_comp = extract_results_by_compartment(Results_extended, dict_comp)


### MASS BALANCE PER COMPARTMENT###

# Estimate mass flows due to the different particle fate process (transfer between compartments, elimination and transformation processes)

# Estimate outflows in mass (g/s) amd number/second
(tables_outputFlows, tables_outputFlows_number) = estimate_outFlows(
    system_particle_object_list, dict_comp
)

# Estimate input flows from transport from other compartments
(tables_inputFlows, tables_inputFlows_num) = estimate_inFlows(
    tables_outputFlows, tables_outputFlows_number, dict_comp, surfComp_list
)

# Create flow dictionaries

# Decode index in input and output flow tables
flows_dict_mass = generate_flows_dict(
    tables_outputFlows, tables_inputFlows, size_dict, MP_form_dict_reverse
)

flows_dict_num = generate_flows_dict(
    tables_outputFlows_number, tables_inputFlows_num, size_dict, MP_form_dict_reverse
)

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


comp_mass_balance_df["Mass balance"] = [
    comp_mass_balance_df["Inflow"][c] - comp_mass_balance_df["Outflow"][c]
    for c in comp_mass_balance_df.index
]


# Add total steady state mass and number of particles concentrations to dataframe


comp_mass_balance_df["Concentration (g/m3)"] = [
    sum(Results_comp_dict[c].concentration_g_m3) for c in comp_mass_balance_df.index
]
comp_mass_balance_df["Concentration (N/m3)"] = [
    sum(Results_comp_dict[c].concentration_num_m3) for c in comp_mass_balance_df.index
]

print(comp_mass_balance_df["Mass balance"])

""" Generate mass and number distribution heatmaps"""


fig_mass, titlename_figmass = plot_fractionDistribution_heatmap(
    Results_extended, fraction="mass_fraction"
)
fig_num, titlename_fignum = plot_fractionDistribution_heatmap(
    Results_extended, fraction="number_fraction"
)


""" Add iput and output flows dict to results extended dataframe"""

Results_extended = addFlows_to_results_df(
    Results_extended, flows_dict_mass, flows_dict_num
)


# This is all in mass units
interactions_pp_df = fillInteractions_fun_OOP_dict(
    system_particle_object_list, SpeciesList, surfComp_list
)

# Create a dictionary of recieving inflows per particle taking the values from the interactions matrix
particle_inflows_dict_mass = {}
particle_inflows_dict_number = {}
for p in system_particle_object_list:
    inflows_p_mass = []
    # inflows_p_num=[]
    for p2 in system_particle_object_list:
        interaction_rate = interactions_pp_df[p2.Pcode][p.Pcode]
        if type(interaction_rate) == dict:
            inflow = {k: v * p2.Pmass_g_SS for k, v in interaction_rate.items()}
            inflows_p_mass.append(inflow)
            # inflows_p_num.append({k: v * p2.Pnum_g_SS for k, v in interaction_rate.items()})
        else:
            inflows_p_mass.append(interaction_rate)
            # inflows_p_num.append(interaction_rate)
    dict_list = [item for item in inflows_p_mass if isinstance(item, dict)]
    # dict_list_num=[item for item in inflows_p_num if isinstance(item, dict)]
    merged_dict = {}
    # merged_dict_num={}
    for d in dict_list:
        for k, v in d.items():
            if k in merged_dict:
                merged_dict[k] += v
                # merged_dict_num[k] += v
            else:
                merged_dict[k] = v
                # merged_dict_num[k] = v

    particle_inflows_dict_mass[p.Pcode] = merged_dict
    # particle_inflows_dict_number[p.Pcode]=merged_dict_num

# Substitute the inputflow values in the results_extended dataframe:

for ele in particle_inflows_dict_mass:
    Results_extended.at[ele, "inflows_g_s"] = particle_inflows_dict_mass[ele]
    # Results_extended.at[ele, "inflows_num_s"] = particle_inflows_dict_number[ele]

Results_extended["Total_inflows_g_s"] = [
    sum(Results_extended.iloc[i].inflows_g_s.values())
    for i in range(len(Results_extended))
]

Results_extended["Total_outflows_g_s"] = [
    sum(Results_extended.iloc[i].outflows_g_s.values())
    for i in range(len(Results_extended))
]

Results_extended["Total_inflows_num_s"] = [
    sum(Results_extended.iloc[i].inflows_num_s.values())
    for i in range(len(Results_extended))
]

Results_extended["Total_outflows_num_s"] = [
    sum(Results_extended.iloc[i].outflows_num_s.values())
    for i in range(len(Results_extended))
]
""" Add iput and output flows dict to compartment results dataframe (results_by_comp)"""
results_by_comp = addFlows_to_results_df_comp(
    results_by_comp, flows_dict_mass, flows_dict_num
)

## Mass and particle number distribution by size fraction
size_distr = [0.5, 5, 50, 500, 5000]
Pmass = []
Pnumber = []
for size in size_distr:
    Pmass.append(
        round(
            sum(
                Results_extended[Results_extended["Size_Fraction_um"] == size]["mass_g"]
            )
            / sum(Results_extended["mass_g"])
            * 100,
            2,
        )
    )
    Pnumber.append(
        round(
            sum(
                Results_extended[Results_extended["Size_Fraction_um"] == size][
                    "number_of_particles"
                ]
            )
            / sum(Results_extended["number_of_particles"])
            * 100,
            2,
        )
    )

# This is the data for overall % in table
size_distribution_df = pd.DataFrame(
    {
        "size_fraction_um": size_distr,
        "percent_of_total_mass": Pmass,
        "percent_of_total_number": Pnumber,
    }
)

""" Estimate exposure indicators """

# For estimating exposure indicators we need to make emissions to targeted compartments.

# Run model with emissions to specific compartments to estimate the emission fractions
from functions.model_run_by_comp import *

model_results = {}
dispersing_comp_list = [
    "Air",
    "Ocean_Mixed_Water",
    "Ocean_Surface_Water",
]

for dispersing_comp in dispersing_comp_list:
    model_results[dispersing_comp] = run_model_comp(
        dispersing_comp,
        input_flow_g_s,
        interactions_df,
        MP_form,
        size_bin,
        particle_forms_coding,
        particle_compartmentCoding,
        system_particle_object_list,
        comp_dict_inverse,
        dict_comp,
        size_dict,
        MP_form_dict_reverse,
        surfComp_list,
    )


#### EXPOSURE INDICATORS ####

# Estimate emission fractions for the chosen emission scenario

emission_fractions_mass_data = emission_fractions_calculations(
    Results_extended,
    model_results,
    dispersing_comp_list,
    dict_comp,
    input_flow_g_s,
    q_num_s,
    size_dict,
    emiss_comp,
)


emiss_fract_fig = plot_emission_fractions(emission_fractions_mass_data, emiss_comp)


# Overall persistance (Pov) and Overall residence time (Tov) in years:
print_output = "True"

(
    Pov_mass_years,
    Pov_num_years,
    Pov_size_dict_years,
    Tov_mass_years,
    Tov_num_years,
    Tov_size_dict_years,
    Pov_Tov_comp_df,
) = Exposure_indicators_calculation(
    tables_outputFlows,
    tables_outputFlows_number,
    Results_extended,
    size_dict,
    dict_comp,
    system_particle_object_list,
    print_output,
    imput_flows_g_s,
    tables_inputFlows_num,
)

""" Add persistence and residence time to results extended dataframe"""

Results_extended = calculate_persistence_residence_time(Results_extended)

# Results by compartment

results_by_comp["Persistence for mass (years)"] = Pov_Tov_comp_df["Pov_years(mass)"]
results_by_comp["Persistence for particle number (years)"] = Pov_Tov_comp_df[
    "Pov_years(particle_number)"
]
results_by_comp["Residence time for mass (years)"] = Pov_Tov_comp_df[
    "Tov_years(mass_g)"
]
results_by_comp["Residence time for particle number (years)"] = Pov_Tov_comp_df[
    "Tov_years(particle_number)"
]


# Caracteristic travel distance (CDT) (m):

# To calculate CTD we need to estimate it by emitting into the especific mobile compartment. We will calculate CTD derived from emmiting to each compartment and taking the higest value:
CTD_mass_list = []
CTD_number_list = []


for CDT_comp in [
    "Ocean_Surface_Water",
    "Ocean_Mixed_Water",
    "Coast_Surface_Water",
    "Coast_Column_Water",
    "Surface_Freshwater",
    "Bulk_Freshwater",
    "Air",
]:
    # input flow (in g per second) for each compartment the User should specify here the input flows per compartment
    q_mass_g_s_dict_CTD = {
        "Ocean_Surface_Water": 0,
        "Ocean_Mixed_Water": 0,
        "Ocean_Column_Water": 0,
        "Coast_Surface_Water": 0,
        "Coast_Column_Water": 0,
        "Surface_Freshwater": 0,
        "Bulk_Freshwater": 0,
        "Sediment_Freshwater": 0,
        "Sediment_Ocean": 0,
        "Sediment_Coast": 0,
        "Beaches_Soil_Surface": 0,
        "Beaches_Deep_Soil": 0,
        "Background_Soil_Surface": 0,
        "Background_Soil": 0,
        "Impacted_Soil_Surface": 0,
        "Impacted_Soil": 0,
        "Air": 0,
    }
    q_mass_g_s_dict_CTD[CDT_comp] = input_flow_g_s

    sp_imputs_CTD = []
    q_mass_g_s_CTD = []
    for compartment in q_mass_g_s_dict_CTD.keys():
        sp_imputs_CTD.append(
            size_bin
            + particle_forms_coding[MP_form]
            + str(particle_compartmentCoding[compartment])
            + "_"
            + boxName
        )
        q_mass_g_s_CTD.append(q_mass_g_s_dict_CTD[compartment])

    imput_flows_g_s_CTD = dict(zip(sp_imputs_CTD, q_mass_g_s_CTD))

    CTD_km = model_run_CTD(
        system_particle_object_list,
        CDT_comp,
        imput_flows_g_s_CTD,
        interactions_df,
        q_mass_g_s_CTD,
        size_dict,
        MP_form_dict_reverse,
        comp_dict_inverse,
        dict_comp,
    )

    CTD_mass_list.append(CTD_km[0])
    CTD_number_list.append(CTD_km[1])

CTD_df = pd.DataFrame(
    index=[
        "Ocean_Surface_Water",
        "Ocean_Mixed_Water",
        "Coast_Surface_Water",
        "Coast_Column_Water",
        "Surface_Freshwater",
        "Bulk_Freshwater",
        "Air",
    ]
)

CTD_df["CTD_mass_km"] = CTD_mass_list
CTD_df["CTD_particle_number_km"] = CTD_number_list

print(
    "Characteristic mass travel distance (CDT): ",
    round(CTD_df["CTD_mass_km"].max(), 1),
    " km",
)

print(
    "Characteristic particle number travel distance (CDT): ",
    round(CTD_df["CTD_particle_number_km"].max(), 1),
    " km",
)

""" Extract input and output flows per compartment """

add_output_flow_conexions(
    results_by_comp,
    dict_comp,
    outputflow_type="outflows_g_s",
    inputflow_type="inflows_g_s",
)
add_output_flow_conexions(
    results_by_comp,
    dict_comp,
    outputflow_type="outflows_num_s",
    inputflow_type="inflows_num_s",
)

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
    results_by_comp,
    tables_outputFlows,
    tables_inputFlows,
    MP_form_dict_reverse,
    size_dict,
    comp_mass_balance_df,
    fig_mass,
    titlename_figmass,
    fig_num,
    titlename_fignum,
    emiss_fract_fig,
)


""" Run Monte Carlo simulation for Sensitivity and Uncertainty Analysis """

# Import Monaco package that runs the Monte Carlo simulations
import monaco as mc

# Import the statistical distributions from scipy.stats that you will be using.
# These must be rv_discrete or rv_continuous functions.
# See https://docs.scipy.org/doc/scipy/reference/stats.html for a complete list.
from scipy.stats import randint, rv_discrete, lognorm, uniform

# Continue code here

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
