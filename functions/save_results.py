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

    # Plot and save RC heatmaps per compartment (Has to be fixed for the procesess where there are multiple values of rate constants)
    # for comp in dict_comp.keys():
    #     plot_heatmap_RC(comp, df4, path_RC)

    # Alternatively we can save the violin plot of rate constant values
    # rateConstants_df=df4
    # selected_columns = rateConstants_df.columns[3:]
    # data_raw=df4[selected_columns]
    # selected_data=data_raw.applymap(sum_if_list)
    # log_data = selected_data.applymap(lambda x: np.log10(x) if x > 0 else np.nan)
    # # Violin Plot
    # plt.figure(figsize=(10, 6))
    # sns.violinplot(data=log_data)
    # #plt.yscale('log')
    # plt.xticks(rotation=90)
    # plt.title("Distribution of rate constants as log(k_s-1)")

    # Save results extended organised in excel sheets per compartment

    results_extended_by_compartment_to_csv(
        path=path_run, results_dict=Results_comp_dict
    )
    ## Plot heatmaps of mass and number distribution

    fig_mass.savefig(path_run + "/" + titlename_figmass + ".png")
    fig_num.savefig(path_run + "/" + titlename_fignum + ".png")

    ## Save Emission Fractions figure

    emiss_fract_fig.savefig(path_run + "/Emission_Fractions.png")

    # Save total mass and number distribution by compartment

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
