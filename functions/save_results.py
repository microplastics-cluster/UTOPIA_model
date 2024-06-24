from datetime import datetime
import os
import matplotlib.pyplot as plt
import seaborn as sns
from functions.plot_results import *


def store_results(
    path,
    outputs_path,
    runName,
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
    fig_mass,
    titlename_figmass,
    fig_num,
    titlename_fignum,
    emiss_fract_fig,
):

    # check whether directory already exists
    if not os.path.exists(path):
        os.mkdir(path)
        print("Folder %s created!" % path)
    else:
        print("Folder %s already exists" % path)

    # check whether directory already exists
    if not os.path.exists(path_run):
        os.mkdir(path_run)
        print("Folder %s created!" % path_run)
    else:
        print("Folder %s already exists" % path_run)

    # plot_rate_constants(df4)

    # # Plot heatmaps of rate constants per compartment

    # Create rate constants folder:
    path_RC = os.path.join(path_run, "RateConstants")

    # check whether directory already exists
    if not os.path.exists(path_RC):
        os.mkdir(path_RC)
        print("Folder %s created!" % path_RC)
    else:
        print("Folder %s already exists" % path_RC)

    # Save rate constants dataframe

    t_filename = os.path.join(path_RC, "RateConstants_table.csv")
    df4.to_csv(t_filename, index=False)

    # Plot and save RC heatmaps per compartment (Has to be fixed for the procesess wherer there are multiple values of rate constants)
    # for comp in dict_comp.keys():
    #     plot_heatmap_RC(comp, df4, path_RC)

    # Save results extended organised in excel sheets per compartment

    results_extended_by_compartment_to_csv(
        path=path_run, results_dict=Results_comp_dict
    )
    ## Plot heatmaps of mass and number distribution

    fig_mass.savefig(path_run + "/" + titlename_figmass + ".png")
    fig_num.savefig(path_run + "/" + titlename_fignum + ".png")

    ## Save Emission Fractions figure

    emiss_fract_fig.savefig(path_run + "/Emission_Fractions.png")

    # Plot results in Total number of particles and Total mass

    # Create Steady State results folders:
    # path_SteadyState_mass = os.path.join(path_run, "SteadyState_mass_distribution")
    # path_SteadyState_number = os.path.join(path_run, "SteadyState_number_distribution")

    # check whether directory already exists
    # if not os.path.exists(path_SteadyState_mass):
    #     os.mkdir(path_SteadyState_mass)
    #     print("Folder %s created!" % path_SteadyState_mass)
    # else:
    #     print("Folder %s already exists" % path_SteadyState_mass)

    # if not os.path.exists(path_SteadyState_number):
    #     os.mkdir(path_SteadyState_number)
    #     print("Folder %s created!" % path_SteadyState_number)
    # else:
    #     print("Folder %s already exists" % path_SteadyState_number)

    # for comp in Results_comp_organiced:
    #     plot_bySize_total_number_particles(
    #         Results_comp_organiced,
    #         comp,
    #         model_lists["dict_size_coding"],
    #         path=path_SteadyState_number,
    #     )
    #     plot_bySize_total_mass(
    #         results_dict=Results_comp_organiced,
    #         comp_name=comp,
    #         dict_size_coding=model_lists["dict_size_coding"],
    #         path=path_SteadyState_mass,
    #     )
    # plot_by(
    #     results_dict=Results_comp_organiced,
    #     comp_name=comp,
    #     dict_size_coding=model_lists["dict_size_coding"],
    #     plot_by="concentration_g_m3",
    # )

    # Save the table of mass distribution
    # massDitribution_filename = os.path.join(
    #     path_SteadyState_mass, "SS_mass_distribution.csv"
    # )
    # df_massDistribution.to_csv(massDitribution_filename)

    # # Save the table of number distribution
    # numberDitribution_filename = os.path.join(
    #     path_SteadyState_number, "SS_number_distribution.csv"
    # )
    # df_numberDistribution.to_csv(numberDitribution_filename)

    # # Plot Mass distribution by compartment and save compartment mass distribution table

    # sns.barplot(data=mass_dist_comp, x="Compartments", y="%_mass").set(
    #     title="Mass distribution"
    # )
    # plt.ylabel("% of total mass")
    # plt.yscale("log")
    # # plt.ylim(0.01, 100)
    # plt.xticks(rotation=90)
    # # Save the plot
    # massDistribPlot_filename = os.path.join(
    #     path_SteadyState_mass, "mass_by_compartment.png"
    # )
    # plt.savefig(massDistribPlot_filename, bbox_inches="tight")
    # plt.show()
    # plt.close()

    # # Plot Number distribution by compartment
    # sns.barplot(data=mass_dist_comp, x="Compartments", y="%_number").set(
    #     title="Particle number distribution"
    # )
    # plt.ylabel("% of total number")
    # plt.yscale("log")
    # # plt.ylim(0.01, 100)
    # plt.xticks(rotation=90)
    # # Save the plot
    # numberDistribPlot_filename = os.path.join(
    #     path_SteadyState_number, "number_by_compartment.png"
    # )
    # plt.savefig(numberDistribPlot_filename, bbox_inches="tight")
    # plt.show()
    # plt.close()

    table_total_mass_number_distribution = os.path.join(
        path_run, "total_distribution_byCompartment.csv"
    )
    mass_dist_comp.to_csv(table_total_mass_number_distribution)

    # Create folder for saving ouput and input mass flows
    path_mass_flows = os.path.join(path_run, "Compartment_mass_flows")

    # check whether directory already exists
    if not os.path.exists(path_mass_flows):
        os.mkdir(path_mass_flows)
        print("Folder %s created!" % path_mass_flows)
    else:
        print("Folder %s already exists" % path_mass_flows)

    for comp in tables_outputFlows:
        flows_tables_to_csv(
            comp,
            tables_outputFlows,
            tables_inputFlows,
            MP_form_dict_reverse,
            size_dict,
            path=path_mass_flows,
        )

    # Save compartment mass balance table
    comp_mass_balance_df.to_csv(os.path.join(path_run, "compartment_mass_balance.csv"))
