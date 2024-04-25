from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.generate_compartmentFlows_tables import *
from helpers.helpers import *
import math


def emission_fractions_calculations(
    Results_extended,
    model_results,
    dispersing_comp_list,
    dict_comp,
    input_flow_g_s,
    q_num_s,
):

    ## Following the LRTP metrics of the emission fractions approach (EFA; φ1, φ2, φ3) from https://doi.org/10.1021/acs.est.2c03047

    ## We use the same values of Crossectional area for Air and Water as in Breivik et al. 2022 and scale it for our water compartments

    Air_crossectional_area_m2 = (
        2.27e9  # Assuming a higth of air of 6000 m (From the OECD tool)
    )
    Water_crossectional_area_m2 = 2.68e7  # Assuming a higth of water of 100 m

    # Assuming that all the water is ocean water.

    Ocean_surface_crosssectional_area_m2 = Water_crossectional_area_m2 * (0.1 / 100)

    Ocean_mixed_crosssectional_area_m2 = Water_crossectional_area_m2 * (
        (100 - 0.1) / 100
    )

    crossSectional_area_m2 = {
        "Air": Air_crossectional_area_m2,
        "Ocean_Surface_Water": Ocean_surface_crosssectional_area_m2,
        "Ocean_Mixed_Water": Ocean_mixed_crosssectional_area_m2,
    }

    """Environmentally Dispersed Fraction (φ1)"""
    # Environmentally Dispersed Fraction (φ1) quantifies the relative extent to which the pollutants (MPs) can reach remote regions.

    φ1_dict_mass = {}
    φ1_dict_num = {}

    # Environmentally Dispersed Fractions (ϕ1)
    for R_comp in dispersing_comp_list:
        Nadv = (
            sum(
                Results_extended["concentration_g_m3"][
                    Results_extended["Compartment"] == R_comp
                ]
            )
            * crossSectional_area_m2[R_comp]
            * float(dict_comp[R_comp].flowVelocity_m_s)
        )
        NE_g_s = input_flow_g_s

        φ1_dict_mass[R_comp] = Nadv / NE_g_s

        Nadv_num = (
            sum(
                Results_extended["concentration_num_m3"][
                    Results_extended["Compartment"] == R_comp
                ]
            )
            * crossSectional_area_m2[R_comp]
            * float(dict_comp[R_comp].flowVelocity_m_s)
        )
        NE_num_s = sum(q_num_s)

        φ1_dict_num[R_comp] = Nadv_num / NE_num_s

    for E1_comp, E1 in zip(φ1_dict_mass.keys(), φ1_dict_mass.values()):
        print(
            "Environmentally Dispersed Mass Fractions through {} = {}".format(
                E1_comp, round(E1, 5)
            )
        )
    print("φ1 for mass =", sum(φ1_dict_mass.values()))

    for E1_comp_num, E1_num in zip(φ1_dict_num.keys(), φ1_dict_num.values()):
        print(
            "Environmentally Dispersed Particle Number Fractions through {} = {}".format(
                E1_comp_num, int(E1_num)
            )
        )
    print("φ1 for particle number =", int(sum(φ1_dict_num.values())))

    """Remotely transferred fraction of mass (ϕ2)"""

    # φ2 expresses the relative extent to which a the MPs are (net) transferred to the target remote compartment following environmental dispersion to the remote region

    internal_comp_process_list = [
        "k_discorporation",
        "k_fragmentation",
        "k_heteroaggregation",
        "k_heteroaggregate_breackup",
        "k_biofouling",
        "k_defouling",
    ]

    φ2_dict_mass = {}
    φ2_dict_num = {}

    # Remotely transferred fraction of mass to the target remote "surface" compartment (φ2) will come through air and water from the compartments listed in the dispersing_comp_list: Air, Coast Surface Water, Coast_Column_Water, Ocean Surface Water and Ocean Mixed Water, :

    """The remotely transferred fraction of particles can only be estimated by size fraction?? If we do it with total number of particles independent of the size it can occur that we get negative values as the flows would be dominated by a specific size fraction and canpotentially occur that the outflows of particles from a compartment in particle number is bigger than the inflows due to the fragmentation processess??"""

    target_remote_comp_List = [
        "Ocean_Surface_Water",
        "Ocean_Column_Water",
        "Sediment_Ocean",
    ]  # ,"Background_Soil_Surface"]

    # target_remote_comp = "Ocean_Surface_Water"
    for target_remote_comp in target_remote_comp_List:
        φ2_Tcomp = {}
        φ2_Tcomp_num = {}
        for transfComp in dispersing_comp_list:
            if transfComp == target_remote_comp:
                NE_x_g_s = NE_g_s
                NE_x_num_s = NE_num_s

            else:
                NE_x_g_s = 0
                NE_x_num_s = 0

            input_flows = model_results[transfComp]["tables_inputFlows"][
                target_remote_comp
            ]
            input_flows_num = model_results[transfComp]["tables_inputFlows_num"][
                target_remote_comp
            ]
            output_flows = model_results[transfComp]["tables_outputFlows"][
                target_remote_comp
            ]
            output_flows_num = model_results[transfComp]["tables_outputFlows_number"][
                target_remote_comp
            ]

            for k_p in internal_comp_process_list:
                if k_p in output_flows:
                    output_flows.drop(columns=k_p, inplace=True)
                if k_p in output_flows_num:
                    output_flows_num.drop(columns=k_p, inplace=True)

            φ2_Tcomp[transfComp] = (
                φ1_dict_mass[transfComp]
                * (
                    NE_x_g_s
                    + sum([sum(input_flows[P]) for P in input_flows])
                    - sum([sum(output_flows[P]) for P in output_flows])
                )
                / NE_g_s
            )
            φ2_Tcomp_num[transfComp] = (
                φ1_dict_num[transfComp]
                * (
                    NE_x_num_s
                    + sum([sum(input_flows_num[P]) for P in input_flows_num])
                    - sum([sum(output_flows_num[P]) for P in output_flows_num])
                )
                / NE_num_s
            )
        φ2_dict_mass[target_remote_comp] = φ2_Tcomp
        φ2_dict_num[target_remote_comp] = φ2_Tcomp_num

    for E2_comp, E2 in zip(φ2_dict_mass.keys(), φ2_dict_mass.values()):
        print(
            "Remotely transferred fraction to {} = {}".format(
                E2_comp, round(sum(E2.values()), 5)
            )
        )

    for E2_compN, E2N in zip(φ2_dict_num.keys(), φ2_dict_num.values()):
        print(
            "Remotely transferred particle number fraction to {} = {}".format(
                E2_compN, sum(E2N.values())
            )
        )

    """Do I need to also consider transfer from ocean column water by making emissions into this compartment?, in this case ϕ1 would be 0 as we consider the flow velocity of the ocean column water to be 0"""

    ## Remotely Accumulated Fraction(φ3): accounts for the degradative loss in surface media and therefore assessess the fraction of the emissions of MPs that accumulates in the target remote surface compartment.

    # Expresses the fraction of deposited chemical that is retained in the respective medium (soil or water) but transferred to deeper layers
    # φ3_Tcomp = {}
    # φ3 = {}
    # for T_comp in φ2_dict["Sediment_Ocean"].keys():

    #     φ3[T_comp] = φ2_dict["Sediment_Ocean"][T_comp] * (
    #         sum(model_results[T_comp]["tables_outputFlows"]["Sediment_Ocean"].k_burial)
    #         / sum(
    #             model_results[T_comp]["tables_outputFlows"][
    #                 "Sediment_Ocean"
    #             ].k_discorporation
    #         )
    #     )

    # φ3_Tcomp["Sediment_Ocean"] = φ3

    # for E3_comp, E3 in zip(φ3_Tcomp.keys(), φ3_Tcomp.values()):
    #     print(
    #         "Remotely Accumulated Fraction to {} = {}".format(E3_comp, sum(E3.values()))
    #     )
