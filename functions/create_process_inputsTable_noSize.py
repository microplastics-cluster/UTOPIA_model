import pandas as pd
import itertools
import os


def write_t_frag_gen_table(model_lists, inputs_path, surfComp_list):
    compNames = model_lists["compartmentNames_list"]
    mpFormsLabels = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
    system_dict = {
        "Compartment": compNames,
        "MPform": mpFormsLabels,
    }

    # Generate all possible combinations
    keys, values = zip(*system_dict.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # Generate dataframe of permutations with input parameter columns
    listOfinputs = [
        "tfrag_gen_d",
    ]

    t_frag_gen_df = pd.DataFrame(permutations_dicts)
    for i in listOfinputs:
        t_frag_gen_df[i] = "NAN"

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
        (t_frag_gen_df["Compartment"] == "Ocean_Surface_Water")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Surface_Water")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Surface_Freshwater")
        & (t_frag_gen_df["MPform"] == "freeMP")
    )
    t_frag_gen_df.loc[cond_frag, "tfrag_gen_d"] = t_frag_gen_FreeSurfaceWater

    cond_frag1 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Surface_Water")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Surface_Water")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Surface_Freshwater")
        & (t_frag_gen_df["MPform"] == "biofMP")
    )
    t_frag_gen_df.loc[cond_frag1, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_biofilm
    )

    cond_frag_new = (
        (t_frag_gen_df["Compartment"] == "Ocean_Surface_Water")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Surface_Water")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Surface_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterMP")
    )

    t_frag_gen_df.loc[cond_frag_new, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_heter
    )

    cond_frag_new1 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Surface_Water")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Surface_Water")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Surface_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
    )

    t_frag_gen_df.loc[cond_frag_new1, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_biofilm * factor_heter
    )

    cond_frag2 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Mixed_Water")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Ocean_Column_Water")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Column_Water")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Bulk_Freshwater")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "freeMP")
    )

    t_frag_gen_df.loc[cond_frag2, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface
    )

    cond_frag3 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Mixed_Water")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Ocean_Column_Water")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Column_Water")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Bulk_Freshwater")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "biofMP")
    )
    t_frag_gen_df.loc[cond_frag3, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface * factor_biofilm
    )

    cond_frag_new2 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Mixed_Water")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Ocean_Column_Water")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Column_Water")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Bulk_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterMP")
    )
    t_frag_gen_df.loc[cond_frag_new2, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_deepWater_soilSurface * factor_heter
    )

    cond_frag_new3 = (
        (t_frag_gen_df["Compartment"] == "Ocean_Mixed_Water")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Ocean_Column_Water")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Coast_Column_Water")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Bulk_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil_Surface")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
    )
    t_frag_gen_df.loc[cond_frag_new3, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater
        * factor_deepWater_soilSurface
        * factor_heter
        * factor_biofilm
    )

    cond_frag4 = (
        (t_frag_gen_df["Compartment"] == "Sediment_Freshwater")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Ocean")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Coast")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil")
        & (t_frag_gen_df["MPform"] == "freeMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil")
        & (t_frag_gen_df["MPform"] == "freeMP")
    )
    t_frag_gen_df.loc[cond_frag4, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment
    )

    cond_frag5 = (
        (t_frag_gen_df["Compartment"] == "Sediment_Freshwater")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Ocean")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Coast")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil")
        & (t_frag_gen_df["MPform"] == "biofMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil")
        & (t_frag_gen_df["MPform"] == "biofMP")
    )
    t_frag_gen_df.loc[cond_frag5, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_biofilm
    )

    cond_frag_new4 = (
        (t_frag_gen_df["Compartment"] == "Sediment_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Ocean")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Coast")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil")
        & (t_frag_gen_df["MPform"] == "heterMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil")
        & (t_frag_gen_df["MPform"] == "heterMP")
    )
    t_frag_gen_df.loc[cond_frag_new4, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_heter
    )

    cond_frag_new5 = (
        (t_frag_gen_df["Compartment"] == "Sediment_Freshwater")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Ocean")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Sediment_Coast")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Urban_Soil")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Background_Soil")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
        | (t_frag_gen_df["Compartment"] == "Agricultural_Soil")
        & (t_frag_gen_df["MPform"] == "heterBiofMP")
    )
    t_frag_gen_df.loc[cond_frag_new5, "tfrag_gen_d"] = (
        t_frag_gen_FreeSurfaceWater * factor_sediment * factor_heter * factor_biofilm
    )

    t_frag_gen_fileName = os.path.join(inputs_path, "t_frag_gen.csv")
    t_frag_gen_df.to_csv(t_frag_gen_fileName)
