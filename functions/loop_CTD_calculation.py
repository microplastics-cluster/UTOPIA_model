from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.generate_compartmentFlows_tables import *
from functions.exposure_indicators_calculation import *


def model_run_CTD(
    system_particle_object_list,
    CDT_comp,
    imput_flows_g_s_CTD,
    interactions_df,
    q_mass_g_s_CTD,
    size_dict,
    MP_form_dict_reverse,
    comp_dict_inverse,
    dict_comp,
):
    """SOLVE SYSTEM OF ODES"""

    R, PartMass_t0 = solve_ODES_SS(
        system_particle_object_list=system_particle_object_list,
        q_num_s=0,
        imput_flows_g_s=imput_flows_g_s_CTD,
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

    Results_extended, mf_shorted, nf_shorted = estimate_fractions(Results)

    # Estimate outflows in mass (g/s) and number/second
    (tables_outputFlows, tables_outputFlows_number) = estimate_outFlows(
        system_particle_object_list, dict_comp
    )

    # Exposure indicators

    # Overall mass persistence

    discorporation_flows = []
    for k in tables_outputFlows:
        discorporation_flows.append(sum(tables_outputFlows[k].k_discorporation))

    Pov_mass_years = (
        (sum(Results_extended["mass_g"]) / sum(discorporation_flows)) / 86400 / 365
    )

    discorporation_flows_num = []
    for k in tables_outputFlows_number:
        discorporation_flows_num.append(
            sum(tables_outputFlows_number[k].k_discorporation)
        )
    Pov_num_years = (
        (sum(Results_extended["number_of_particles"]) / sum(discorporation_flows_num))
        / 86400
        / 365
    )

    # Caracteristic travel distance

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

    return (CTD_mass_m / 1000, CTD_number_m / 1000)
