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
from functions.save_results import *
from functions.generate_compartmentFlows_tables import *
from functions.exposure_indicators_calculation import *


def model_run(
    model_lists,
    comp_dict_inverse,
    particle_compartmentCoding,
    system_particle_object_list,
    interactions_df,
    q_mass_g_s,
    imput_flows_g_s,
    particle_forms_coding,
    size_dict,
    MP_form_dict_reverse,
    dict_comp,
    surfComp_list,
    compartment,
    saveName,
    saveOpt,
    df4,
):
    """SOLVE SYSTEM OF ODES"""

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

    # print("Distribution of mass in the system")
    # print(mf_shorted[:10])
    df_massDistribution = mf_shorted[:10]

    # print("distribution of particle number in the system")
    # print(nf_shorted[:10])
    df_numberDistribution = nf_shorted[:10]

    # Mass distribution by compartment
    mass_frac_100 = []
    num_frac_100 = []
    mass_conc_g_m3 = []
    num_conc = []
    for comp in list(dict_comp.keys()):
        mass_frac_100.append(
            sum(
                Results_extended[Results_extended["Compartment"] == comp][
                    "mass_fraction"
                ]
            )
            * 100
        )
        num_frac_100.append(
            sum(
                Results_extended[Results_extended["Compartment"] == comp][
                    "number_fraction"
                ]
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

    # Estimate outflows in mass (g/s) and number/second
    (tables_outputFlows, tables_outputFlows_number) = estimate_outFlows(
        system_particle_object_list, dict_comp
    )

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
    # print(comp_mass_balance_df)

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
        sum(Results_comp_dict[c].concentration_num_m3)
        for c in comp_mass_balance_df.index
    ]

    """ Estimate exposure indicators """

    # Exposure indicators
    (
        Pov_mass_years,
        Pov_num_years,
        Pov_size_dict_sec,
        Tov_mass_years,
        Tov_num_years,
        Tov_size_dict_sec,
    ) = Exposure_indicators_calculation(
        tables_outputFlows,
        tables_outputFlows_number,
        Results_extended,
        size_dict,
        dict_comp,
    )
    # Caracteristic travel distance

    CTD_mass_km, CTD_number_km = calculate_CTD(
        Pov_mass_years, Results_extended, dict_comp, Pov_num_years, compartment
    )

    """ Generate mass and number distribution heatmaps"""

    plot_fractionDistribution_heatmap(Results_extended, fraction="mass_fraction")
    plot_fractionDistribution_heatmap(Results_extended, fraction="number_fraction")

    # Save results

    if saveOpt == "save":

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
    else:
        pass
