from functions.solver_SteadyState import *
from functions.extract_results import *
from functions.generate_compartmentFlows_tables import *
from helpers.helpers import *
import math


def emission_fractions_calculations(
    Results_extended,
    model_results_size_bin,
    dispersing_comp_list,
    dict_comp,
    input_flow_g_s,
    q_num_s,
    size_bin,
    size_dict,
    emiss_comp,
):
    model_results = model_results_size_bin[size_bin]
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

    """ Mass Emission Fractions"""

    """Environmentally Dispersed Fraction (φ1)"""
    # Environmentally Dispersed Fraction (φ1) quantifies the relative extent to which the pollutants (MPs) can reach remote regions.

    φ1_dict_mass = {}

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

    # Estimate composition of Environmentally Dispersed Fractions in percentage:
    φ1_mass_comp = [
        round(v * 100 / sum(φ1_dict_mass.values()), 2) for v in φ1_dict_mass.values()
    ]

    for E1_comp, E1 in zip(φ1_dict_mass.keys(), φ1_dict_mass.values()):
        print(
            "Environmentally Dispersed Mass Fractions through {} = {}".format(
                E1_comp, round(E1, 5)
            )
        )
    print("φ1 for mass =", round(sum(φ1_dict_mass.values()), 5))

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

    # Remotely transferred fraction of mass to the target remote compartment (φ2) will come through air and water from the compartments listed in the dispersing_comp_list: Air, Ocean_Mixed_Water and Ocean_Surface_Water.

    # The target remote compartments are Ocean Surface Water as an approximation to study transfer to the Ocean Gyres, Ocean Column water and Ocean sediment and background soil surface for representing transer to remote beaches.

    target_remote_comp_List = [
        "Ocean_Surface_Water",
        "Ocean_Column_Water",
        "Sediment_Ocean",
        "Background_Soil_Surface",
    ]

    for target_remote_comp in target_remote_comp_List:
        φ2_Tcomp = {}
        for transfComp in dispersing_comp_list:
            if transfComp == target_remote_comp:
                NE_x_g_s = NE_g_s
            else:
                NE_x_g_s = 0

            input_flows = model_results[transfComp]["tables_inputFlows"][
                target_remote_comp
            ]

            output_flows = model_results[transfComp]["tables_outputFlows"][
                target_remote_comp
            ]

            for k_p in internal_comp_process_list:
                if k_p in output_flows:
                    output_flows.drop(columns=k_p, inplace=True)

            # Substitute the columns of the dataframe that have list of values for the sum of the values(sum output flows of that process)
            for k_o in output_flows:
                output_flows[k_o] = [
                    sum(x) if isinstance(x, list) else x for x in output_flows[k_o]
                ]

            φ2_Tcomp[transfComp] = (
                φ1_dict_mass[transfComp]
                * (
                    NE_x_g_s
                    + sum([sum(input_flows[P]) for P in input_flows])
                    - sum([sum(output_flows[P]) for P in output_flows])
                )
                / NE_g_s
            )

        φ2_dict_mass[target_remote_comp] = φ2_Tcomp

    φ2_mass = []
    for E2_comp, E2 in zip(φ2_dict_mass.keys(), φ2_dict_mass.values()):
        φ2_mass.append(sum(E2.values()))

        print(
            "Remotely transferred fraction to {} = {}".format(
                E2_comp, round(sum(E2.values()), 5)
            )
        )

    print(
        "Total remotely transferred mass fraction = {}".format(round(sum(φ2_mass), 5))
    )

    φ2_mass_table = pd.DataFrame(
        {"Remotely transferred fraction to": φ2_dict_mass.keys(), "φ2": φ2_mass}
    )

    emission_fractions_mass_data = {
        "Emission Fraction": ["φ1", "φ2_1", "φ2_2", "φ2_3", "φ2_4"],
        "y": [round(sum(φ1_dict_mass.values()), 5)] + φ2_mass,
    }

    plot_emission_fractions(emission_fractions_mass_data, emiss_comp)

    """ Particle Number Emission Fractions are estimated per size fraction"""

    """Environmentally Dispersed Fraction (φ1)"""
    # Environmentally Dispersed Fraction (φ1) quantifies the relative extent to which the pollutants (MPs) can reach remote regions.

    φ1_dict_num = {}
    for size, s_num in zip(size_dict.keys(), size_dict.values()):
        φ1_dict_num[size] = {}
        for R_comp in dispersing_comp_list:
            Nadv_num = (
                sum(
                    Results_extended["concentration_num_m3"][
                        Results_extended["Compartment"] == R_comp
                    ][Results_extended["Size_Fraction_um"] == s_num]
                )
                * crossSectional_area_m2[R_comp]
                * float(dict_comp[R_comp].flowVelocity_m_s)
            )
            NE_num_s = sum(q_num_s)

            φ1_dict_num[size][R_comp] = Nadv_num / NE_num_s

        # for E1_comp_num, E1_num in zip(φ1_dict_num[size].keys(), φ1_dict_num[size].values()):
        #     print(
        #     "Environmentally Dispersed fraction of "+str(s_num)+" µm MP particles through {} = {}".format(
        #         E1_comp_num, int(E1_num)
        #     )
        # )
        print(
            "φ1 for particle number for " + str(s_num) + " µm MP particles =",
            int(sum(φ1_dict_num[size].values())),
        )

        print(pd.DataFrame(φ1_dict_num))

        # φ1_num_comp = [
        #     round(v * 100 / sum(φ1_dict_num.values()), 2) for v in φ1_dict_num.values()
        # ]

        # for E1_comp_num, E1_num in zip(φ1_dict_num.keys(), φ1_dict_num.values()):
        #     print(
        #         "Environmentally Dispersed Particle Number Fractions through {} = {}".format(
        #             E1_comp_num, int(E1_num)
        #         )
        #     )
        # print("φ1 for particle number =", int(sum(φ1_dict_num.values())))

        """Remotely transferred fraction of particle number (ϕ2)"""

        # φ2 expresses the relative extent to which a the MPs are (net) transferred to the target remote compartment following environmental dispersion to the remote region

        # We can estimate φ2 in mass idependent of the size fraction of the particles, however when estimating φ2 in particle number we have to do it per size fraction.
        """The remotely transferred fraction of particle number can only be estimated by size fraction?? If we do it with total number of particles independent of the size it can occur that we get negative values as the flows would be dominated by a specific size fraction and can potentially occur that the outflows of particles from a compartment in particle number is bigger than the inflows due to the fragmentation processess??"""

        φ2_dict_num = {}
        for size, s_num in zip(size_dict.keys(), size_dict.values()):
            φ2_dict_num[size] = {}
            for target_remote_comp in target_remote_comp_List:
                φ2_Tcomp_num = {}
                for transfComp in dispersing_comp_list:
                    if transfComp == target_remote_comp:
                        NE_x_num_s = NE_num_s
                    else:
                        NE_x_num_s = 0

                    input_flows_num = model_results_size_bin[size][transfComp][
                        "tables_inputFlows_num"
                    ][target_remote_comp]

                    output_flows_num = model_results_size_bin[size][transfComp][
                        "tables_outputFlows_number"
                    ][target_remote_comp]

                    for k_p in internal_comp_process_list:
                        if k_p in output_flows_num:
                            output_flows_num.drop(columns=k_p, inplace=True)

                    # Substitute the columns of the dataframe that have list of values for the sum of the values(sum output flows of that process)
                    for k_o in output_flows_num:
                        output_flows_num[k_o] = [
                            sum(x) if isinstance(x, list) else x
                            for x in output_flows_num[k_o]
                        ]

                    φ2_Tcomp_num[transfComp] = (
                        φ1_dict_num[size][transfComp]
                        * (
                            NE_x_num_s
                            + sum([sum(input_flows_num[P]) for P in input_flows_num])
                            - sum([sum(output_flows_num[P]) for P in output_flows_num])
                        )
                        / NE_num_s
                    )
                φ2_dict_num[size][target_remote_comp] = φ2_Tcomp_num

            φ2_mass = []
            for E2_comp, E2 in zip(φ2_dict_mass.keys(), φ2_dict_mass.values()):
                φ2_mass.append(sum(E2.values()))

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

            print(
                "Total remotely transferred mass fraction = {}".format(
                    round(sum(φ2_mass), 5)
                )
            )

            φ2_mass_table = pd.DataFrame(
                {"Remotely transferred fraction to": φ2_dict_mass.keys(), "φ2": φ2_mass}
            )


def plot_emission_fractions(emission_fractions_data, emiss_comp):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    emission_fractions_data["y"] = [np.log10(x) for x in emission_fractions_data["y"]]

    data = {
        "x": ["φ1", "φ2_1", "φ2_2", "φ2_3", "φ2_4"],
        "y": emission_fractions_data["y"],
        "category": [
            "φ1",
            "φ2_Ocean_Surface_Water",
            "φ2_Ocean_Column_Water",
            "φ2_Sediment_Ocean",
            "φ2_Background_Soil_Surface",
        ],
    }
    df = pd.DataFrame(data)

    # Create a dictionary to map unique categories to colors
    colors = {
        "φ1": "red",
        "φ2_Ocean_Surface_Water": "blue",
        "φ2_Ocean_Column_Water": "green",
        "φ2_Sediment_Ocean": "orange",
        "φ2_Background_Soil_Surface": "purple",
    }

    # Plot the DataFrame with different colors per category and add a legend
    fig, ax = plt.subplots()
    for category, color in colors.items():
        df_category = df[df["category"] == category]
        ax.scatter(
            df_category["x"], df_category["y"], c=color, label=category, marker="_"
        )

    plt.ylabel("log10(φ)")

    plt.title("Mass Emission Fractions after emissions to " + emiss_comp)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.show()

    # emission_fractions_df=pd.DataFrame(emission_fractions_data)
    # # Plot the DataFrame with logarithm of y-axis values
    # emission_fractions_df.plot( x="Emission Fraction", y="y",kind='scatter',marker='_')
    # #plt.xlabel('Remote Compartment')
    # plt.ylabel('log10(φ)')

    # plt.show()
    pass
