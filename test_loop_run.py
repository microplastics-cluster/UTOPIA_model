import os
from datetime import datetime
import pandas as pd
import numpy as np
from model_run import *


inputs_path = os.path.join(os.path.dirname(__file__), "inputs")


"""Define run parameters"""

## Define microplastics physical properties

# The user can also select a preloaded file instead of typing in the values. In this case the user wont need to run the code between lines 29 and 34 and neither the code between lines 42 and 50. The user will have to run line 56 with the selected input file


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
mp_imputFile_name = os.path.join(
    inputs_path, "inputs_microplastics.csv"
)  # Choose one existing input file to load

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
q_mass_g_s = 1

# particle imput
size_bin = "e"
MP_form = "freeMP"
MP_density = "lowDensity"  # To be changed based on the MP imputs file

for comp in compList:
    compartment = comp

    saveName = (
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
        particle_forms_coding,
        size_dict,
        MP_form_dict_reverse,
        size_bin,
        compartment,
        MP_form,
        saveName,
    )