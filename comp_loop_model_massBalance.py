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


def model_run(
    mp_imputFile_name,
    comp_impFile_name,
    comp_interactFile_name,
    q_mass_g_s,
    size_bin,
    compartment,
    MP_form,
):

    inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

    """Define run parameters"""

    if mp_imputFile_name == "\inputs_microplastics.csv":
        MP_density = "lowDensity"
    elif mp_imputFile_name == "\inputs_microplastics_HDP.csv":
        MP_density = "highDensity"
    else:
        print("define scenario name based on particle density choosing a MP_density")

    boxName = "Utopia"

    runName = (
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

    """Generate objects"""

    # Generate objects
    (
        modelBoxes,
        system_particle_object_list,
        SpeciesList,
        process_inputs_df,
        spm,
        dict_comp,
        model_lists,
        particles_df,
        MPforms_list,
    ) = generate_objects(
        inputs_path,
        boxName=boxName,
        comp_impFile_name=comp_impFile_name,
        comp_interactFile_name=comp_interactFile_name,
        mp_imputFile_name=mp_imputFile_name,
    )

    surfComp_list = [c for c in dict_comp if "Surface" in c]

    """Estimate rate constants per particle"""

    for particle in system_particle_object_list:
        generate_rateConstants(particle, spm, dict_comp)

    #### Original results with no time limit

    ## create rate constants table:
    RC_df = create_rateConstants_table(system_particle_object_list)
    df4 = RC_df.fillna(0)

    # Plot rate constants

    """(FIX RC for wet deposition, now its given as a list of rate constants per surface compartment only for dry deposition)This needs to be fixed also for the matrix of interactions and estimation of flows"""

    """Build Matrix of interactions"""

    interactions_df = fillInteractions_fun_OOP(
        system_particle_object_list, SpeciesList, surfComp_list
    )

    from functions.fill_interactions_Knames import *

    interactions_df_Knames = fillInteractions_Knames(
        system_particle_object_list, SpeciesList
    )

    """SOLVE SYSTEM OF ODES"""

    particle_compartmentCoding = dict(
        zip(
            model_lists["compartmentNames_list"],
            list(range(len(model_lists["compartmentNames_list"]))),
        )
    )
    comp_dict_inverse = {v: k for k, v in particle_compartmentCoding.items()}

    sp_imput = (
        size_bin
        + particle_forms_coding[MP_form]
        + str(particle_compartmentCoding[compartment])
        + "_"
        + boxName
    )

    print(
        "Imput flow = "
        + str(q_mass_g_s)
        + " g/s of "
        + MP_form
        + " of size "
        + str(size_dict[size_bin])
        + " into the "
        + compartment
        + " compartment"
    )

    R, PartMass_t0 = solve_ODES_SS(
        system_particle_object_list=system_particle_object_list,
        q_mass_g_s=q_mass_g_s,
        q_num_s=0,
        sp_imput=sp_imput,
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

    mass_dist_comp = pd.DataFrame(columns=["Compartments"])
    mass_dist_comp["Compartments"] = list(dict_comp.keys())
    mass_dist_comp["%_mass"] = mass_frac_100
    mass_dist_comp["%_number"] = num_frac_100

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
        sum(Results_comp_dict[c].concentration_num_m3)
        for c in comp_mass_balance_df.index
    ]

    """ Estimate exposure indicators """

    # Overall residence time
    overall_residence_time_calculation(tables_outputFlows, Results_extended)
