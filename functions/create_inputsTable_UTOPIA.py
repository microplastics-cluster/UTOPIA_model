# creates a pandas dataframe of process parameters inputs for all the particles in the system
# (regarding combination of sizes, MPforms and compartments)

# If the systems consist on several model boxes a further
# element has to be added to the function

import pandas as pd
import itertools


def create_inputsTable_UTOPIA(
    inputs_path, model_lists, thalf_deg_d_dict, alpha_hetr_dict
):
    compNames = model_lists["compartmentNames_list"]
    mpFormsLabels = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
    # sizeBins = ["x01um", "um", "x10um", "x100um", "mm"]
    sizeBinsLables = list(model_lists["dict_size_coding"].keys())

    system_dict = {
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
    # Assumptions:
    # Heteroaggregated particles degrade 10 times slower than the free MPs
    # Biofouled particles degrade 5 times slower than the free MPs

    for key in thalf_deg_d_dict:
        cond = dataFrame_inputs["MPform"] == key
        dataFrame_inputs.loc[cond, "thalf_deg_d"] = thalf_deg_d_dict[key]

    # Timescale for fragmentation of the 1000um size fraction (mp5): tfrag_gen_d

    "Old Assumption (Full Multi): fragmentation only occurs for free and biofouled MPs and the timescale depends on the compartment and aggregation state"
    "In UTOPIA we include fragmentation of the heteroaggregated MPs as being 100 slower than fragmentation of the Free MPs and breackup of biofouled and heteroaggregated will be two times slowed of those only heteroaggregated, following the same assumption as for free and biofouled. These values are used in the Domercq et al. 2021 paper and they are asumptions made from lack of current knowlegde"  #!Values to be revisited

    ## Assumptions:

    # fractionation does not take place in the Air compartment. -->To be revisited!!
    # fragmentation is fastest when the particles are in free form and in the surface water compartments :fragmentation of Free MPs in the surface water compartments takes 36.5 days to occur
    # fragemntation of biofouled particles takes double the time than for Free particles and for heteroaggregated particles it takes 100 times more
    # Fragmentation in the lower water compartments and in the surface of the sediment takes 10 times more time than in the surface water
    # Fragmentation in the sediment compartments take 100 times more time than in the surface water compartments

    t_frag_gen_FreeSurfaceWater = 36.5
    factor_biofilm = 2
    factor_heter = 100
    factor_deepWater_soilSurface = 10
    factor_sediment = 100

    cond_frag = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        & (dataFrame_inputs["MPform"] == "freeMP")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        & (dataFrame_inputs["MPform"] == "freeMP")
    )

    dataFrame_inputs.loc[cond_frag, "tfrag_gen_d"] = t_frag_gen_FreeSurfaceWater

    cond_frag1 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag1, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_biofilm
    )

    cond_frag_new = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )

    dataFrame_inputs.loc[cond_frag_new, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_heter
    )

    cond_frag_new1 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Surface_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Surface_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )

    dataFrame_inputs.loc[cond_frag_new1, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_biofilm * factor_heter
    )

    cond_frag2 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )

    dataFrame_inputs.loc[cond_frag2, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface
    )

    cond_frag3 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag3, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface * factor_biofilm
    )

    cond_frag_new2 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag_new2, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface * factor_heter
    )

    cond_frag_new3 = (
        (dataFrame_inputs["Compartment"] == "Ocean_Mixed_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Ocean_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Coast_Column_Water")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Bulk_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil_Surface")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag_new3, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater
        * factor_deepWater_soilSurface
        * factor_heter
        * factor_biofilm
    )

    cond_frag4 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "freeMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag4, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment
    )

    cond_frag5 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "biofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag5, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_biofilm
    )

    cond_frag_new4 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "heterMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag_new4, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_heter
    )

    cond_frag_new5 = (
        (dataFrame_inputs["Compartment"] == "Sediment_Freshwater")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Ocean")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Sediment_Coast")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Urban_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Background_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
        | (dataFrame_inputs["Compartment"] == "Agricultural_Soil")
        & (dataFrame_inputs["MPform"] == "heterBiofMP")
        & (dataFrame_inputs["sizeBin"] == "mp5")
    )
    dataFrame_inputs.loc[cond_frag_new5, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_heter * factor_biofilm
    )

    # Time for the biofilm coverage to grow: tbiof_growth_d

    "here we follow the hypothesis that biofouling occurs at slower rates in deeper waters due to reduced light limiting the growth of the biofilm organisms (Kooi et al., 2017).Values of time for biofim growth are based on experimental findings that indicate that biofilm formation takes place within days or weeks (Rummel et al., 2017)."

    # To be implemented: Product Formulation Controls the Impact of Biofouling on Consumer Plastic Photochemical Fate in the Ocean (Nelson et al. 2021)

    # Biofouling is modelled to occur in free and heteroaggregated particles

    # Assumptions
    # biofilm growth is slow in deeper waters. Biofilm growth takes 10 days in surface water compartments and 30 days in the ocean mixed waters and coast column water and 300 days in the deep ocean

    tbiof_growth_surfaceWater_d = 10
    tbiof_growth_lowDeepWater_d = 30
    tbiof_growth_deepWater_d = 300

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

    dataFrame_inputs.loc[cond_biof1, "tbiof_growth_d"] = tbiof_growth_surfaceWater_d

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
    dataFrame_inputs.loc[cond_biof2, "tbiof_growth_d"] = tbiof_growth_lowDeepWater_d

    cond_biof3 = (dataFrame_inputs["MPform"] == "freeMP") & (
        dataFrame_inputs["Compartment"] == "Ocean_Column_Water"
    ) | (dataFrame_inputs["MPform"] == "heterMP") & (
        dataFrame_inputs["Compartment"] == "Ocean_Column_Water"
    )

    dataFrame_inputs.loc[cond_biof3, "tbiof_growth_d"] = tbiof_growth_deepWater_d

    # Defouling (and its time rate measure tbiof_degrade_d) is the disintegration of the biofilm layer.

    "it can occur due to light limitation, grazing, or dissolution of carbonates in acid waters (Kooi et al., 2017).So far assumed as null due to lack of data regarding biofilm degradation times."

    # Defouling would be only modelled for the biofouled particles (biofMP and heterBiofMP?) To be decided if its depth dependent also (therefore compartment dependent)

    # Heteroaggregation attachment efficiency: alpha_heter.

    "Heteroaggegation happens to free and biofouled particles. It is hypothesized that biofilm increases the attachment efficiency of a plastic particle, reflected in two times higher values of  for biofiouled plastic particles compared to the pristine form. We assumed there is no heteroaggregation in the sediment or any soil compartment and neither in air"
    # REF value: Besseling et al. 2017

    alpha_heter_Free = float(alpha_hetr_dict["freeMP"])
    alpha_heter_biof = float(alpha_hetr_dict["biofMP"])

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
    dataFrame_inputs.loc[cond_alpha1, "alpha_heter"] = alpha_heter_Free

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
    dataFrame_inputs.loc[cond_alpha2, "alpha_heter"] = alpha_heter_biof

    # Output dataFrame_inputs as csv file

    dataFrame_inputs.to_csv(inputs_path + "\processInputs_table.csv", index=False)

    return dataFrame_inputs


"""List of references for the parameterization of the environmental compartments characteristics"""

# def create_compartment_inputsTable():

#     compList = [
#         "Ocean_Surface_Water",
#         "Ocean_Mixed_Water",
#         "Ocean_Column_Water",
#         "Coast_Surface_Water",
#         "Coast_Column_Water",
#         "Surface_Freshwater",
#         "Bulk_Freshwater",
#         "Sediment_Freshwater",
#         "Sediment_Ocean",
#         "Sediment_Coast",
#         "Urban_Soil_Surface",
#         "Urban_Soil",
#         "Background_Soil_Surface",
#         "Background_Soil",
#         "Agricultural_Soil_Surface",
#         "Agricultural_Soil",
#         "Air",
#     ]

#     earth_SA_m2 = 5.10e14
#     land_SA_m2 = earth_SA_m2 * 0.29
#     water_SA_m2 = earth_SA_m2 * 0.71
#     """Total Surface Area of the earth and water and land percentages taken form OECD tool"""

#     freshWater_SA_m2 = water_SA_m2 * 0.025
#     oceanSeaWater_SA_m2 = water_SA_m2 * 0.975
#     """"""

#     coastWater_SA_m2 = oceanSeaWater_SA_m2 * 0.025
#     oceanWater_SA_m2 = oceanSeaWater_SA_m2 * 0.975

#     agri_land_SA_m2 = (
#         land_SA_m2 * 0.38
#     )  # Ref: FAO 2020. Land use in agriculture by the numbers (https://www.fao.org/sustainability/news/detail/en/c/1274219/)
#     urban_land_SA_m2 = (
#         land_SA_m2 * 0.02
#     )  # Liu, Z., He, C., Zhou, Y. and Wu, J., 2014. How much of the world’s land has been urbanized, really? A hierarchical framework for avoiding confusion. Landscape Ecology, 29, pp.763-771.
#     background_land_SA_m2 = land_SA_m2 * 0.60

#     flow_velocity_ocean_m_s = 0.002  # Ref: from The OECD Pov and LRTP Screening Tool (Version 2.2). F. Wegmann et al(2009), Environmental Modeling & Software 24, 228-237.

#     freshWater_discharge_km3_yr = 37288  # The river-based estimate of global continental discharge presented here is 37 288 ± 662 km3 yr−1. Ref: Dai, A. and Trenberth, K. (2002). Estimates of Freshwater Discharge from Continents: Latitudinal and Seasonal Variations. Journal of Hydrometeorology 3(6) pp. 660-687.

#     depth_coastSurface_water_m = 0.1
#     continental_shell_depth_m = 50
#     continental_shell_width_m = 60000

#     waterFlow_surface_freshWater_m3_s = (
#         freshWater_discharge_km3_yr * 1e9 / (60 * 60 * 24 * 365)
#     )
#     waterFlow_coastColumn_water_m3_s = (
#         (coastWater_SA_m2 / continental_shell_width_m) * continental_shell_depth_m
#     ) * flow_velocity_ocean_m_s
#     waterFlow_coastSurface_water_m3_s = (
#         (coastWater_SA_m2 / continental_shell_width_m) * depth_coastSurface_water_m
#     ) * flow_velocity_ocean_m_s
#     waterFlow_oceanSurface_water_m3_s = waterFlow_coastSurface_water_m3_s
#     waterFlow_oceanMixed_water_m3_s = waterFlow_coastColumn_water_m3_s
