# creates a pandas dataframe of process parameters inputs for all the particles in the system
# (regarding combination of sizes, MPforms and compartments)

# If the systems consist on several model boxes a further
# element has to be added to the function

import pandas as pd
import itertools


def create_inputsTable_UTOPIA(compartments, modelBoxes, inputs_path):
    boxNames = [b.Bname for b in modelBoxes]
    compNames = [c.Cname for c in compartments]
    mpFormsLabels = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
    # sizeBins = ["x01um", "um", "x10um", "x100um", "mm"]
    sizeBinsLables = ["mp5", "mp4", "mp3", "mp2", "mp1"]

    # system_dict = {
    #     "Compartment": compNames,
    #     "MPform": mpFormsLabels,
    #     "sizeBin": sizeBinsLables,
    # }
    system_dict = {
        "modelBox": boxNames,
        "Compartment": compNames,
        "MPform": mpFormsLabels,
        "sizeBin": sizeBinsLables,
    }

    # Generate all possible combinations
    keys, values = zip(*system_dict.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # Generate dataframe of permutations with input parameter columns
    listOfinputs = [
        "thalf_deg_d",
        "tfrag_gen_d",
        "tbiof_growth_d",
        "tbiof_degrade_d",
        "alpha_heter",
    ]
    dataFrame_inputs = pd.DataFrame(permutations_dicts)
    for i in listOfinputs:
        dataFrame_inputs[i] = "NAN"

    # Stablish input parameter values

    ## Degradation half time: thalf_deg_d
    "Values used in Domercq et al. 2021, go to publication for more details on the selection of these values and asumptions made"

    thalf_deg_d_dict = {
        "freeMP": 5000,
        "heterMP": 50000,
        "biofMP": 25000,
        "heterBiofMP": 100000,
    }

    for key in thalf_deg_d_dict:
        cond = dataFrame_inputs["MPform"] == key
        dataFrame_inputs.loc[cond, "thalf_deg_d"] = thalf_deg_d_dict[key]

    # Timescale for fragmentation of the 1000um size fraction (mp1): tfrag_gen_d

    "Old Assumption (Full Multi): fragmentation only occurs for free and biofouled MPs and the timescale depends on the compartment and aggregation state"
    "In UTOPIA we include fragmentation of the heteroaggregated MPs as being 100 slower than fragmentation of the Free MPs and breackup of biofouled and heteroaggregated will be two times slowed of tose only heteroaggregated, following the same assumption as for free and biofouled. These values are used in the Domercq et al. 2021 paper and they are asumptions made from lack of current knowlegde"  #!Values to be revisited

    # Assumptions: fractionation does not take place in the Air compartment. -->To be revisited!!

    cond_frag = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        & (dataFrame_inputs["MPform"] == "freeMP")
    )

    dataFrame_inputs.loc[cond_frag, "tfrag_gen_d"] = 36.5

    cond_frag1 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag1, "tfrag_gen_d"] = 73

    cond_frag_new = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )

    dataFrame_inputs.loc[cond_frag_new, "tfrag_gen_d"] = 36.5 * 100

    cond_frag_new1 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )

    dataFrame_inputs.loc[cond_frag_new1, "tfrag_gen_d"] = 36.5 * 100 * 2

    cond_frag2 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )

    dataFrame_inputs.loc[cond_frag2, "tfrag_gen_d"] = 365

    cond_frag3 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag3, "tfrag_gen_d"] = 730

    cond_frag_new2 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag_new2, "tfrag_gen_d"] = 365 * 100

    cond_frag_new3 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag_new3, "tfrag_gen_d"] = 365 * 100 * 2

    cond_frag4 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag4, "tfrag_gen_d"] = 3650

    cond_frag5 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag5, "tfrag_gen_d"] = 7300

    cond_frag_new4 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag_new4, "tfrag_gen_d"] = 3650 * 100

    cond_frag_new5 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp1")
    )
    dataFrame_inputs.loc[cond_frag_new5, "tfrag_gen_d"] = 3650 * 100 * 2

    # Time for the biofilm coverage to grow: tbiof_growth_d

    "here we follow the hypothesis that biofouling occurs at slower rates in deeper waters due to reduced light limiting the growth of the biofilm organisms (Kooi et al., 2017).Values of time for biofim growth are based on experimental findings that indicate that biofilm formation takes place within days or weeks (Rummel et al., 2017)."

    # To be implemented: Product Formulation Controls the Impact of Biofouling on Consumer Plastic Photochemical Fate in the Ocean (Nelson et al. 2021)

    # Biofouling is modelled to occur in free and heteroaggregated particles#

    cond_biof1 = (
        (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        | (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
    )

    dataFrame_inputs.loc[cond_biof1, "tbiof_growth_d"] = 10

    cond_biof2 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
    )
    dataFrame_inputs.loc[cond_biof2, "tbiof_growth_d"] = 30

    cond_biof3 = (dataFrame_inputs["MPform"] == "freeMP") & (
        dataFrame_inputs["Compartment"] == "Ocean_Column_Water"
    ) | (dataFrame_inputs["MPform"] == "heterMP") & (
        dataFrame_inputs["Compartment"] == "Ocean_Column_Water"
    )

    dataFrame_inputs.loc[cond_biof3, "tbiof_growth_d"] = 300

    # Defouling (and its time rate measure tbiof_degrade_d) is the disintegration of the biofilm layer.

    "it can occur due to light limitation, grazing, or dissolution of carbonates in acid waters (Kooi et al., 2017).So far assumed as null due to lack of data regarding biofilm degradation times."

    # Defouling would be only modelled for the biofouled particles (biofMP and heterBiofMP?) To be decided if its depth dependent also (therefore compartment dependent)

    # Heteroaggregation attachment efficiency: alpha_heter.

    "Heteroaggegation happens to free and biofouled particles. It is hypothesized that biofilm increases the attachment efficiency of a plastic particle, reflected in two times higher values of  for biofiouled plastic particles compared to the pristine form. We assumed there is no heteroaggregation in the sediment or any soil compartment and neither in air"
    # REF value: Besseling et al. 2017

    cond_alpha1 = (dataFrame_inputs["MPform"] == "freeMP") & (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        | (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
    )
    dataFrame_inputs.loc[cond_alpha1, "alpha_heter"] = 0.01

    cond_alpha2 = (dataFrame_inputs["MPform"] == "biofMP") & (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        | (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
    )
    dataFrame_inputs.loc[cond_alpha2, "alpha_heter"] = 0.02

    # Output dataFrame_inputs as csv file

    dataFrame_inputs.to_csv(inputs_path + "\processInputs_table.csv", index=False)

    return dataFrame_inputs


def create_compartment_inputsTable():

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


earth_SA_m2 = 5.10e14
land_SA_m2 = earth_SA_m2 * 0.29
water_SA_m2 = earth_SA_m2 * 0.71
"""Total Surface Area of the earth and water and land percentages taken form OECD tool"""

freshWater_SA_m2 = water_SA_m2 * 0.025
oceanSeaWater_SA_m2 = water_SA_m2 * 0.975
""""""

coastWater_SA_m2 = oceanSeaWater_SA_m2 * 0.025
oceanWater_SA_m2 = oceanSeaWater_SA_m2 * 0.975

agri_land_SA_m2 = land_SA_m2 * 0.02
urban_land_SA_m2 = land_SA_m2 * 0.38
background_land_SA_m2 = land_SA_m2 * 0.60
