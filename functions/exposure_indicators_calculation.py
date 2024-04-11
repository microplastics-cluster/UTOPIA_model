import pandas as pd

size_list = ["a", "b", "c", "d", "e"]


def Exposure_indicators_calculation(
    tables_outputFlows,
    tables_outputFlows_number,
    Results_extended,
    size_dict,
    dict_comp,
):
    """Overall persistance (years)"""

    # Overall persistance for the plastic material in all size classes:

    # Overall mass persistence

    discorporation_flows = []
    for k in tables_outputFlows:
        discorporation_flows.append(sum(tables_outputFlows[k].k_discorporation))

    Pov_mass_sec = sum(Results_extended["mass_g"]) / sum(discorporation_flows)

    Pov_mass_days = Pov_mass_sec / 86400
    Pov_mass_years = Pov_mass_days / 365

    print("Overall mass persistence (years): " + str(Pov_mass_years))

    # Overall number persistence

    discorporation_flows_num = []
    for k in tables_outputFlows_number:
        discorporation_flows_num.append(
            sum(tables_outputFlows_number[k].k_discorporation)
        )

    Pov_num_sec = sum(Results_extended["number_of_particles"]) / sum(
        discorporation_flows_num
    )

    Pov_num_days = Pov_num_sec / 86400
    Pov_num_years = Pov_num_days / 365

    print("Overall particle number persistence (years): " + str(Pov_num_years))

    # Overall persistence specific to each size class are mass and number independent:

    size_list = ["a", "b", "c", "d", "e"]
    Pov_size_dict_sec = {}
    for size in size_list:
        discorporation_fargmentation_flows = []
        for k in tables_outputFlows:
            outputFlow_df = tables_outputFlows[k]
            sliced_df = outputFlow_df[outputFlow_df.index.str[0] == size]
            discorporation_fargmentation_flows.append(
                sum(sliced_df.k_fragmentation) + sum(sliced_df.k_discorporation)
            )

        mass_sizeFraction = sum(
            Results_extended[Results_extended.index.str[0] == size].mass_g
        )

        Pov_size_sec = mass_sizeFraction / sum(discorporation_fargmentation_flows)
        Pov_size_days = Pov_size_sec / 86400
        Pov_size_years = Pov_size_days / 365
        print(
            "Overall persistence of size "
            + str(size_dict[size])
            + " um (years): "
            + str(Pov_size_years)
        )
        Pov_size_dict_sec[size] = Pov_size_sec

    """ Overall residence time (years)"""  # Includes burial in the deep sediment and will also inlcude sequestration in the deep soils once the transport process in the soil compartments are included

    systemloss_flows_mass = []
    systemloss_flows_number = []

    for k in tables_outputFlows:
        # if "Soil" in k:
        #     systemloss_flows.append(
        #         sum(tables_outputFlows[k].k_sequestration_deep_soils)
        #     )
        if "Sediment" in k:
            systemloss_flows_mass.append(
                sum(tables_outputFlows[k].k_discorporation)
                + sum(tables_outputFlows[k].k_burial)
            )
            systemloss_flows_number.append(
                sum(tables_outputFlows_number[k].k_discorporation)
                + sum(tables_outputFlows_number[k].k_burial)
            )
        else:
            systemloss_flows_mass.append(sum(tables_outputFlows[k].k_discorporation))
            systemloss_flows_number.append(
                sum(tables_outputFlows_number[k].k_discorporation)
            )

    Tov_mass_sec = sum(Results_extended["mass_g"]) / sum(systemloss_flows_mass)
    Tov_num_sec = sum(Results_extended["number_of_particles"]) / sum(
        systemloss_flows_number
    )

    Tov_mass_days = Tov_mass_sec / 86400
    Tov_num_days = Tov_num_sec / 86400
    Tov_mass_years = Tov_mass_days / 365
    Tov_num_years = Tov_num_days / 365

    print("Overall mass residence time (years): " + str(Tov_mass_years))
    print("Overall particle number residence time (years): " + str(Tov_num_years))

    # Overall residence time specific to each size class (mass and number independent):
    Tov_size_dict_sec = {}
    for size in size_list:
        systemloss_flows_size = []
        for k in tables_outputFlows:
            outputFlow_df = tables_outputFlows[k]
            sliced_df = outputFlow_df[outputFlow_df.index.str[0] == size]
            # if "Soil" in k:
            #     systemloss_flows_size.append(
            #         sum(sliced_df.k_fragmentation)
            #         + sum(sliced_df.k_discorporation)
            #         + sum(sliced_df.k_sequestration_deep_soils)
            #     )
            if "Sediment" in k:
                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation)
                    + sum(sliced_df.k_discorporation)
                    + sum(sliced_df.k_burial)
                )
            else:
                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation) + sum(sliced_df.k_discorporation)
                )

        mass_sizeFraction = sum(
            Results_extended[Results_extended.index.str[0] == size].mass_g
        )

        Tov_size_sec = mass_sizeFraction / sum(systemloss_flows_size)
        Tov_size_days = Tov_size_sec / 86400
        Tov_size_years = Tov_size_days / 365
        print(
            "Overall residence time of size "
            + str(size_dict[size])
            + " um (years): "
            + str(Tov_size_years)
        )
        Tov_size_dict_sec[size] = Tov_size_sec

    return (
        Pov_mass_years,
        Pov_num_years,
        Pov_size_dict_sec,
        Tov_mass_years,
        Tov_num_years,
        Tov_size_dict_sec,
    )


import os
from datetime import datetime
import pandas as pd
import numpy as np
from model_run import *
from functions.generate_MPinputs_table import *


def calculate_CTD(Pov_mass_years, Results_extended, dict_comp, Pov_num_years, CDT_comp):
    """Characteristic travel distance (CDT) (m)"""

    # Characteristic travel distance for for the plastic material in all size classes. We do not calculate CDT for the ocean column water as we do not have flow velocity data on this (most probably has a much slower flow velocity than the surface and mixed ocean water)

    # CDT of particles mass and number

    CTD_mass_m = (
        (Pov_mass_years * 365 * 24 * 60 * 60)
        * sum(Results_extended[Results_extended["Compartment"] == CDT_comp].mass_g)
        / sum(Results_extended["mass_g"])
        * float(dict_comp[CDT_comp].flowVelocity_m_s)
    )
    CTD_number_m = (
        (Pov_num_years * 365 * 24 * 60 * 60)
        * sum(
            Results_extended[
                Results_extended["Compartment"] == CDT_comp
            ].number_of_particles
        )
        / sum(Results_extended["number_of_particles"])
        * float(dict_comp[CDT_comp].flowVelocity_m_s)
    )

    print(
        "Characteristic travel distance (CDT) of particles mass (km) for "
        + CDT_comp
        + ": "
        + str(CTD_mass_m / 1000)
        + " km"
    )

    print(
        "Characteristic travel distance (CDT) of particles number (km) for "
        + CDT_comp
        + ": "
        + str(CTD_number_m / 1000)
        + " km"
    )

    return (CTD_mass_m / 1000, CTD_number_m / 1000)

    # CTD_mass_dic_km = {}
    # CTD_number_dic_km = {}
    # for m in [
    #     "Ocean_Surface_Water",
    #     "Ocean_Mixed_Water",
    #     "Coast_Surface_Water",
    #     "Coast_Column_Water",
    #     "Surface_Freshwater",
    #     "Bulk_Freshwater",
    #     "Air",
    # ]:
    #     CTD_mass_m = (
    #         (Pov_mass_years * 365 * 24 * 60 * 60)
    #         * sum(Results_extended[Results_extended["Compartment"] == m].mass_g)
    #         / sum(Results_extended["mass_g"])
    #         * float(dict_comp[m].flowVelocity_m_s)
    #     )
    #     CTD_number_m = (
    #         (Pov_num_years * 365 * 24 * 60 * 60)
    #         * sum(
    #             Results_extended[
    #                 Results_extended["Compartment"] == m
    #             ].number_of_particles
    #         )
    #         / sum(Results_extended["number_of_particles"])
    #         * float(dict_comp[m].flowVelocity_m_s)
    #     )
    #     CTD_mass_km = CTD_mass_m / 1000
    #     CTD_mass_dic_km[m] = CTD_mass_km
    #     CTD_number_dic_km[m] = CTD_number_m / 1000

    #     data_CTD = {
    #         "CTD_mass_km": CTD_mass_dic_km.values(),
    #         "CTD_number_km": CTD_number_dic_km.values(),
    #     }
    #     CTD_df = pd.DataFrame(data=data_CTD, index=CTD_mass_dic_km.keys())

    #     print(CTD_df)

    #     # print("Characteristic mass travel distance for " + m + " (km): " + str(CTD_mass_km))
    #     # print("Characteristic particle number travel distance for " + m + " (km): " + str(CTD_number_m/1000))

    # Characteristic travel distance per size class

    # CTD_data_dic_km = {}
    # for m in [
    #     "Ocean_Surface_Water",
    #     "Ocean_Mixed_Water",
    #     "Coast_Surface_Water",
    #     "Coast_Column_Water",
    #     "Surface_Freshwater",
    #     "Bulk_Freshwater",
    #     "Air",
    # ]:
    #     Results_extended_m = Results_extended[Results_extended["Compartment"] == m]
    #     CTD_size_km = []
    #     for size in size_list:
    #         CTD_size_m = (
    #             Pov_size_dict_sec[size]
    #             * sum(
    #                 Results_extended_m[Results_extended_m.index.str[0] == size].mass_g
    #             )
    #             / sum(Results_extended[Results_extended.index.str[0] == size].mass_g)
    #             * float(dict_comp[m].flowVelocity_m_s)
    #         )
    #         CTD_size_km.append(CTD_size_m / 1000)
    #     CTD_data_dic_km[m] = CTD_size_km
    # CTD_size_df_km = pd.DataFrame(
    #     CTD_data_dic_km, index=[size_dict[i] for i in size_list]
    # )
    # print(CTD_size_df_km)
