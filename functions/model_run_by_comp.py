from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.generate_compartmentFlows_tables import *
from helpers.helpers import *
import math


def run_model_comp(
    recieving_comp,
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
):
    """Run the model to solve the system of ODEs when emissions are made to the specified compartment: recieving_comp"""
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
    q_mass_g_s_dict[recieving_comp] = input_flow_g_s

    """SOLVE SYSTEM OF ODES"""

    sp_imputs = []
    q_mass_g = []
    for compartment in q_mass_g_s_dict.keys():
        sp_imputs.append(
            size_bin
            + particle_forms_coding[MP_form]
            + str(particle_compartmentCoding[compartment])
            + "_"
            + "Utopia"
        )
        q_mass_g.append(q_mass_g_s_dict[compartment])

    imput_flows_g_s = dict(zip(sp_imputs, q_mass_g))

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

    # Test that there are no negative results
    for i, idx in zip(R["mass_g"], R.index):
        if i < 0:
            print("negative values in the solution for " + idx)
        else:
            pass

    # Estimate mass and number fractions and extract ranking tables of the species with higest fractions to understand the distribution of the particles in the system by mass and number of particles

    Results_extended, mf_shorted, nf_shorted = estimate_fractions(Results)

    # Estimate mass flows due to the different particle fate process (transfer between compartments, elimination and transformation processes)

    # Estimate outflows in mass (g/s) amd number/second
    (tables_outputFlows, tables_outputFlows_number) = estimate_outFlows(
        system_particle_object_list, dict_comp
    )

    # Estimate imput flows from transport from other compartments
    (tables_inputFlows, tables_inputFlows_num) = estimate_inFlows(
        tables_outputFlows, tables_outputFlows_number, dict_comp, surfComp_list
    )

    return {
        "Results_extended": Results_extended,
        "tables_outputFlows": tables_outputFlows,
        "tables_outputFlows_number": tables_outputFlows_number,
        "tables_inputFlows": tables_inputFlows,
        "tables_inputFlows_num": tables_inputFlows_num,
    }
