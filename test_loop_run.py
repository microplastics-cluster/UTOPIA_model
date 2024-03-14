import os
from datetime import datetime
import pandas as pd
import numpy as np
from model_run import *
from functions.generate_MPinputs_table import *


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


## Suspended particulates properties

spm_diameter_um = 0.5
spm_density_kg_m3 = 2000

## Select fragmentation type

frag_styles_dict = {
    "sequential_fragmentation": np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0.1, 0.9, 0, 0, 0],
            [0.01, 0.09, 0.9, 0, 0],
            [0.001, 0.009, 0.09, 0.9, 0],
        ]
    ),
    "erosive_fragmentation": np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0.9, 0.1, 0, 0, 0],
            [0.99, 0.009, 0.001, 0, 0],
            [0.99, 0.009, 0.0009, 0.0001, 0],
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

frag_style = "erosive_fragmentation"

fsd = frag_styles_dict[frag_style]

# optionally the user can type its own fsd matrix following the desciption below:

# estimate fragmentation relation between size bins using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html)
# each article fractions into fragments of samller sizes and the distribution is expresses via the fragment size distribution matrix fsd. # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class

## choose input files to load

comp_impFile_name = "\inputs_compartments.csv"
comp_interactFile_name = (
    "\compartment_interactions.csv"  # Fixed, should not be modified
)
# mp_imputFile_name = os.path.join(
#     inputs_path, "inputs_microplasticsPE.csv"
# )  # Choose one existing input file to load

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
MPforms_list = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))

MP_form_dict_reverse = {v: k for k, v in particle_forms_coding.items()}

## Compartments:

compList = [
    "Ocean_Surface_Water",
    "Ocean_Mixed_Water",
    "Ocean_Column_Water",
    "Coast_Surface_Water",
    "Coast_Column_Water",
    "Surface_Freshwater",
    "Bulk_Freshwater",
    "Sediment_Freshwater",
    "Sediment_Ocean",
    "Sediment_Coast",
    "Urban_Soil_Surface",
    "Urban_Soil",
    "Background_Soil_Surface",
    "Background_Soil",
    "Agricultural_Soil_Surface",
    "Agricultural_Soil",
    "Air",
]

# input flow (in g per second)


# particle imput
size_bin = "e"
MP_form = "freeMP"
MP_density = "lowDensity"  # To be changed based on the MP imputs file
input_flow_g_s = 1

# save results option
saveOpt = "Notsave"  # "save" or "Notsave"


for comp in compList:
    compartment = comp
    # input flow (in g per second) for each compartment the User should specify here the input flows per compartment
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
        "Urban_Soil_Surface": 0,
        "Urban_Soil": 0,
        "Background_Soil_Surface": 0,
        "Background_Soil": 0,
        "Agricultural_Soil_Surface": 0,
        "Agricultural_Soil": 0,
        "Air": 0,
    }
    q_mass_g_s_dict[compartment] = input_flow_g_s

    print("Input flows for ", compartment, " in g per second: ", input_flow_g_s)

    saveName = (
        MP_density
        + "MP_Emissions_"
        + MP_form
        + "_"
        + str(size_dict[size_bin])
        + "_nm_"
        + compartment
    )

    particle_compartmentCoding = dict(
        zip(
            compList,
            list(range(len(compList))),
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

    model_run(
        inputs_path,
        boxName,
        MPforms_list,
        mp_imputFile_name,
        comp_impFile_name,
        comp_interactFile_name,
        spm_diameter_um,
        spm_density_kg_m3,
        fsd,
        q_mass_g_s,
        imput_flows_g_s,
        particle_forms_coding,
        size_dict,
        MP_form_dict_reverse,
        size_bin,
        compartment,
        MP_form,
        saveName,
        saveOpt,
    )
