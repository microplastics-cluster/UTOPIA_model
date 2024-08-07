import pandas as pd

size_list = ["a", "b", "c", "d", "e"]


def Exposure_indicators_calculation(
    tables_outputFlows,
    tables_outputFlows_number,
    Results_extended,
    size_dict,
    dict_comp,
    system_particle_object_list,
):
    #### EXPOSURE INDICATORS ####
    # When estimating the exposure indicators we do not take on account the Column Water and Ocean Sediment compartments. This is to mantain consistency with the OECD tool as there the particles going deeper than 100 m into the ocean are considered lossess, therefore we also use that as a boundary in our system. Also in this way we prevent the ocean sediment and column water from driving the POV and residence times values. However our emission fraction estimates do take this compartmets into consideration and the MPs fate into the whole UTOPIA systme is reflected there.

    """Overall persistance (years)"""

    # Overall persistance for the plastic material in all size classes and table of overall persistance per compartment:

    # Overall mass persistence

    discorporation_flows = []
    discorporation_flows_all = []
    for k in tables_outputFlows:
        discorporation_flows_all.append(sum(tables_outputFlows[k].k_discorporation))
        if k != "Ocean_Column_Water" and k != "Sediment_Ocean":
            discorporation_flows.append(sum(tables_outputFlows[k].k_discorporation))

    # remove ocean column water and ocean sediment results
    comp_outBoundares = ["Ocean_Column_Water", "Sediment_Ocean"]
    Results_extended_EI = Results_extended[
        ~Results_extended["Compartment"].isin(comp_outBoundares)
    ]

    Pov_mass_sec = sum(Results_extended_EI["mass_g"]) / sum(discorporation_flows)

    Pov_mass_days = Pov_mass_sec / 86400
    Pov_mass_years = Pov_mass_days / 365

    print("Overall mass persistence (years): " + str(int(Pov_mass_years)))

    # Table of overall mass persistences per compartment
    Pov_comp_years = []
    for c in tables_outputFlows.keys():
        if sum(Results_extended[Results_extended["Compartment"] == c]["mass_g"]) == 0:
            Pov_comp_years.append(
                "NaN"
            )  # When there is no mass in the compartment Pov has no value (marked as NAN)
        else:
            Pov_comp_years.append(
                sum(Results_extended[Results_extended["Compartment"] == c]["mass_g"])
                / sum(tables_outputFlows[c].k_discorporation)
                / 86400
                / 365
            )

    Pov_Tov_comp_df = pd.DataFrame(
        {"Compartment": tables_outputFlows.keys(), "Pov_years(mass)": Pov_comp_years}
    )

    # Overall number persistence

    discorporation_flows_num = []
    discorporation_flows_num_all = []
    for k in tables_outputFlows_number:
        discorporation_flows_num_all.append(
            sum(tables_outputFlows_number[k].k_discorporation)
        )
        if k != "Ocean_Column_Water" and k != "Sediment_Ocean":
            discorporation_flows_num.append(
                sum(tables_outputFlows_number[k].k_discorporation)
            )

    Pov_num_sec = sum(Results_extended_EI["number_of_particles"]) / sum(
        discorporation_flows_num
    )

    Pov_num_days = Pov_num_sec / 86400
    Pov_num_years = Pov_num_days / 365

    print("Overall particle number persistence (years): " + str(int(Pov_num_years)))

    # Table of discorporation flows per compartment in number
    Pov_comp_years_num = []
    for c in tables_outputFlows.keys():
        if (
            sum(
                Results_extended[Results_extended["Compartment"] == c][
                    "number_of_particles"
                ]
            )
            == 0
        ):
            Pov_comp_years_num.append(
                "NaN"
            )  # When there is no number in the compartment Pov has no value (marked as NAN)
        else:
            Pov_comp_years_num.append(
                sum(
                    Results_extended[Results_extended["Compartment"] == c][
                        "number_of_particles"
                    ]
                )
                / sum(tables_outputFlows_number[c].k_discorporation)
                / 86400
                / 365
            )

    Pov_Tov_comp_df["Pov_years(particle_number)"] = Pov_comp_years_num

    # Overall persistence specific to each size class are mass and number independent:

    # NOTE! When the mass is only present in one size fraction then the Pov has to be equal to the overall Pov and mas and number Pov should be the same

    size_list = ["a", "b", "c", "d", "e"]
    Pov_size_dict_sec = {}
    for size in size_list:
        discorporation_fargmentation_flows = []
        for k in tables_outputFlows:
            if k != "Ocean_Column_Water" and k != "Sediment_Ocean":
                outputFlow_df = tables_outputFlows[k]
                sliced_df = outputFlow_df[outputFlow_df.index.str[0] == size]
                discorporation_fargmentation_flows.append(
                    sum(sliced_df.k_fragmentation) + sum(sliced_df.k_discorporation)
                )

        mass_sizeFraction = sum(
            Results_extended_EI[Results_extended_EI.index.str[0] == size].mass_g
        )
        if (
            mass_sizeFraction == 0
        ):  ## If there are no particles of a specific size fraction in the system one does not need to estimate Pov
            Pov_size_dict_sec[size] = "NaN"
            continue

        Pov_size_sec = mass_sizeFraction / sum(discorporation_fargmentation_flows)
        Pov_size_days = Pov_size_sec / 86400
        Pov_size_years = Pov_size_days / 365
        print(
            "Overall persistence of size "
            + str(size_dict[size])
            + " um (years): "
            + str(int(Pov_size_years))
        )
        Pov_size_dict_sec[size] = Pov_size_sec

    """ Overall residence time (years)"""
    # With the new system boundaries wew acount for sequestration in deep soils and burial into coast and frehwater sediment but for the Ocean sediment we do not take burial but the settling into the ocean column water compartment as well as mixing. (We exclude the Ocean column water and ocean sediment from the system boundaries in these calculations)

    systemloss_flows_mass = []
    systemloss_flows_number = []

    for k in tables_outputFlows:
        if k in ["Urban_Soil", "Background_Soil", "Agricultural_Soil"]:
            systemloss_flows_mass.append(
                sum(tables_outputFlows[k].k_discorporation)
                + sum(tables_outputFlows[k].k_sequestration_deep_soils)
            )
        if k in ["Sediment_Freshwater", "Sediment_Coast"]:
            systemloss_flows_mass.append(
                sum(tables_outputFlows[k].k_discorporation)
                + sum(tables_outputFlows[k].k_burial)
            )
            systemloss_flows_number.append(
                sum(tables_outputFlows_number[k].k_discorporation)
                + sum(tables_outputFlows_number[k].k_burial)
            )
        elif k == "Ocean_Mixed_Water":
            # Calculate net flow lossess from the Ocean Midex water compartment through deep ocean
            flow_mix_down_mass = []
            flow_mix_down_number = []
            flow_mix_up_mass = []
            flow_mix_up_number = []
            flow_rising_mass = []
            flow_rising_number = []
            for p in system_particle_object_list:
                if p.Pcompartment.Cname == "Ocean_Mixed_Water":
                    flow_mix_down_mass.append(
                        p.RateConstants["k_mixing"][1] * p.Pmass_g_SS
                    )
                    flow_mix_down_number.append(
                        p.RateConstants["k_mixing"][1] * p.Pnum_SS
                    )
                elif p.Pcompartment.Cname == "Ocean_Column_Water":
                    flow_mix_up_mass.append(p.RateConstants["k_mixing"] * p.Pmass_g_SS)
                    flow_mix_up_number.append(p.RateConstants["k_mixing"] * p.Pnum_SS)
                    flow_rising_mass.append(p.RateConstants["k_rising"] * p.Pmass_g_SS)
                    flow_rising_number.append(p.RateConstants["k_rising"] * p.Pnum_SS)

            systemloss_flows_mass.append(
                sum(tables_outputFlows[k].k_discorporation)
                + sum(tables_outputFlows[k].k_settling)
                + sum(flow_mix_down_mass)
                - sum(flow_mix_up_mass)
                - sum(flow_rising_mass)
            )
            systemloss_flows_number.append(
                sum(tables_outputFlows_number[k].k_discorporation)
                + sum(tables_outputFlows_number[k].k_settling)
                + sum(flow_mix_down_number)
                - sum(flow_mix_up_number)
                - sum(flow_rising_number)
            )

        elif k != "Ocean_Column_Water" and k != "Sediment_Ocean":
            systemloss_flows_mass.append(sum(tables_outputFlows[k].k_discorporation))
            systemloss_flows_number.append(
                sum(tables_outputFlows_number[k].k_discorporation)
            )
        else:
            pass

    Tov_mass_sec = sum(Results_extended_EI["mass_g"]) / sum(systemloss_flows_mass)
    Tov_num_sec = sum(Results_extended_EI["number_of_particles"]) / sum(
        systemloss_flows_number
    )

    Tov_mass_days = Tov_mass_sec / 86400
    Tov_num_days = Tov_num_sec / 86400
    Tov_mass_years = Tov_mass_days / 365
    Tov_num_years = Tov_num_days / 365

    print(
        "Overall residence time is calculated assuming the model boundaries to be at 100 m depth into the Ocean, 30 cm into the sediments and 0.1 m into the soil. Particles travelling deeper are considered losses"
    )

    print("Overall mass residence time (years): " + str(round(Tov_mass_years, 1)))

    print(
        "Overall particle number residence time (years): "
        + str(round(Tov_num_years, 1))
    )

    # NOTE: When only one size class pressent, should the residence time be the same in particle number and in mass?? !!! TO check!!!

    # Overall residence time specific to each size class (mass and number independent):

    Tov_size_dict_sec = {}
    for size in size_list:

        mass_sizeFraction = sum(
            Results_extended_EI[Results_extended_EI.index.str[0] == size].mass_g
        )

        if mass_sizeFraction == 0:
            Tov_size_dict_sec[size] = "NaN"
            continue

        systemloss_flows_size = []
        for k in tables_outputFlows:
            outputFlow_df = tables_outputFlows[k]
            sliced_df = outputFlow_df[outputFlow_df.index.str[0] == size]
            if k in ["Urban_Soil", "Background_Soil", "Agricultural_Soil"]:
                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation)
                    + sum(sliced_df.k_discorporation)
                    + sum(sliced_df.k_sequestration_deep_soils)
                )
            elif k in ["Sediment_Freshwater", "Sediment_Coast"]:
                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation)
                    + sum(sliced_df.k_discorporation)
                    + sum(sliced_df.k_burial)
                )
            elif k == "Ocean_Mixed_Water":
                # Calculate net flow lossess from the Ocean Midex water compartment through deep ocean
                flow_mix_down = []
                flow_mix_up = []
                flow_rising = []
                for p in system_particle_object_list:
                    if p.Pcode[0] == size:
                        if p.Pcompartment.Cname == "Ocean_Mixed_Water":
                            flow_mix_down.append(
                                p.RateConstants["k_mixing"][1] * p.Pmass_g_SS
                            )
                        elif p.Pcompartment.Cname == "Ocean_Column_Water":
                            flow_mix_up.append(
                                p.RateConstants["k_mixing"] * p.Pmass_g_SS
                            )
                            flow_rising.append(
                                p.RateConstants["k_rising"] * p.Pmass_g_SS
                            )

                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation)
                    + sum(sliced_df.k_discorporation)
                    + sum(sliced_df.k_settling)
                    + sum(flow_mix_down)
                    - sum(flow_mix_up)
                    - sum(flow_rising)
                )

            elif k != "Ocean_Column_Water" and k != "Sediment_Ocean":
                systemloss_flows_size.append(
                    sum(sliced_df.k_fragmentation) + sum(sliced_df.k_discorporation)
                )

            else:
                pass

        Tov_size_sec = mass_sizeFraction / sum(systemloss_flows_size)
        Tov_size_days = Tov_size_sec / 86400
        Tov_size_years = Tov_size_days / 365
        print(
            "Overall residence time of size "
            + str(size_dict[size])
            + " um (years): "
            + str(round(Tov_size_years, 2))
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
from functions.model_run import *
from functions.generate_MPinputs_table import *


def calculate_CTD(Pov_mass_years, Results_extended, dict_comp, Pov_num_years, CDT_comp):
    """Characteristic travel distance (CDT) (m)"""

    # Characteristic travel distance for for the plastic material in all size classes. We do not calculate CDT for the ocean column water as we do not have flow velocity data on this (most probably has a much slower flow velocity than the surface and mixed ocean water)

    # CDT of particles mass and number

    # remove ocean column water and ocean sediment results (new model boundaries)
    comp_outBoundares = ["Ocean_Column_Water", "Sediment_Ocean"]
    Results_extended_EI = Results_extended[
        ~Results_extended["Compartment"].isin(comp_outBoundares)
    ]

    CTD_mass_m = (
        (Pov_mass_years * 365 * 24 * 60 * 60)
        * sum(
            Results_extended_EI[Results_extended_EI["Compartment"] == CDT_comp].mass_g
        )
        / sum(Results_extended_EI["mass_g"])
        * float(dict_comp[CDT_comp].flowVelocity_m_s)
    )
    CTD_number_m = (
        (Pov_num_years * 365 * 24 * 60 * 60)
        * sum(
            Results_extended_EI[
                Results_extended_EI["Compartment"] == CDT_comp
            ].number_of_particles
        )
        / sum(Results_extended_EI["number_of_particles"])
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
