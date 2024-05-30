import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
import os

import pandas as pd

particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


# Define a function to check if a value is a list and calculate the sum if it is
def sum_if_list(value):
    if isinstance(value, list):
        return sum(value)
    else:
        return value


def plot_bySize_total_number_particles(results_dict, comp_name, dict_size_coding, path):
    new_size_dict = dict(
        zip(
            [particle_sizes_coding[x] for x in dict_size_coding],
            [str(y) for y in dict_size_coding.values()],
        )
    )

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15, 7), sharey=True)

    for i, agg in enumerate(results_dict[comp_name]):
        y = results_dict[comp_name][agg]["number_of_particles"]
        X = results_dict[comp_name][agg]["species"]
        X_1 = [s[0] for s in X]
        x = [new_size_dict[s] for s in X_1]
        L = list(zip([float(z) for z in x], y))
        L.sort()
        L_shorted = [(str(x), y) for (x, y) in L]
        xs = [x for x, y in L_shorted]
        ys = [y for x, y in L_shorted]
        axs[i].bar(xs, ys)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel("Total Number of Particles")
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)

    # Save the plots
    totalNumber_filename = os.path.join(
        path, comp_name + "_SS_totalNumber_distribution.png"
    )
    plt.savefig(totalNumber_filename, bbox_inches="tight")

    plt.close(fig)


def plot_bySize_total_mass(results_dict, comp_name, dict_size_coding, path):
    new_size_dict = dict(
        zip(
            [particle_sizes_coding[x] for x in dict_size_coding],
            [str(y) for y in dict_size_coding.values()],
        )
    )

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15, 7), sharey=True)

    for i, agg in enumerate(results_dict[comp_name]):
        y = results_dict[comp_name][agg]["mass_g"]
        X = results_dict[comp_name][agg]["species"]
        X_1 = [s[0] for s in X]
        x = [new_size_dict[s] for s in X_1]
        L = list(zip([float(z) for z in x], y))
        L.sort()
        L_shorted = [(str(x), y) for (x, y) in L]
        xs = [x for x, y in L_shorted]
        ys = [y for x, y in L_shorted]
        axs[i].bar(xs, ys)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel("Total mass (g)")
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)

    # Save the plots
    totalMass_filename = os.path.join(
        path, comp_name + "_SS_totalMass_distribution.png"
    )
    plt.savefig(totalMass_filename, bbox_inches="tight")

    plt.close(fig)


def plot_by(results_dict, comp_name, dict_size_coding, plot_by):
    new_size_dict = dict(
        zip(
            [particle_sizes_coding[x] for x in dict_size_coding],
            [str(y) for y in dict_size_coding.values()],
        )
    )

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15, 7), sharey=True)

    for i, agg in enumerate(results_dict[comp_name]):
        y = results_dict[comp_name][agg][plot_by]
        X = results_dict[comp_name][agg]["species"]
        X_1 = [s[0] for s in X]
        x = [new_size_dict[s] for s in X_1]
        L = list(zip([float(z) for z in x], y))
        L.sort()
        L_shorted = [(str(x), y) for (x, y) in L]
        xs = [x for x, y in L_shorted]
        ys = [y for x, y in L_shorted]
        axs[i].bar(xs, ys)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel(plot_by)
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)

    plt.close(fig)

    return fig


def extract_output_table(
    comp, tables_outputFlows, MP_form_dict_reverse, size_dict, path
):
    T = tables_outputFlows[comp]
    MP_size = []
    MP_form = []
    for x in T.index:
        MP_size.append(size_dict[x[0]])
        MP_form.append(MP_form_dict_reverse[x[1:2]])
    T.insert(0, "MP_size", MP_size)
    T.insert(1, "MP_form", MP_form)

    # Save table
    outputFlows_filename = os.path.join(path, comp + "_output_mass_flows.csv")
    T.to_csv(outputFlows_filename)

    # # Print heatmaps
    # if sum(T.loc[:, T.columns[2:]].max()) == 0:
    #     print("All values are 0 and no heatmap can be printed")
    # else:
    #     ax = plt.axes()
    #     sns.heatmap(
    #         T.loc[:, T.columns[2:]],
    #         xticklabels=True,
    #         yticklabels=True,
    #         norm=LogNorm(),
    #         linewidths=1,
    #         linecolor="grey",
    #         ax=ax,
    #     )
    #     plt.title("Output Flows for " + comp + " (g/s)")
    #     plt.xlabel("Process")
    #     plt.ylabel("Particle")
    #     plt.show()
    return T


def sum_if_list(value):
    if isinstance(value, list):
        return sum(value)
    else:
        return value


def plot_heatmap_RC(comp, rateConstants_df, path_RC):
    T = rateConstants_df[rateConstants_df["Compartment"] == comp]

    selected_columns = T.columns[4:]
    column1 = T["MP_form"]
    column2 = T["Size_Bin"]

    # Select the desired columns
    df_selected = T[selected_columns]
    df_selected = df_selected.applymap(sum_if_list)

    # Generate the heatmap with combined labels
    fig, ax = plt.subplots(
        figsize=(len(df_selected.columns) * 0.7, len(df_selected) * 0.4)
    )
    sns.heatmap(
        df_selected,
        xticklabels=df_selected.columns,
        yticklabels=column1 + "_" + column2,
        norm=LogNorm(),
        linewidths=1,
        linecolor="grey",
        ax=ax,
    )

    # Save the heatmap as an image
    heatmap_filename = os.path.join(path_RC, comp + "_RC_heatmap.png")
    plt.savefig(heatmap_filename, bbox_inches="tight")

    # Add a title to the heatmap
    ax.set_title(comp + " rate constants (s-1)")

    plt.close(fig)


def flows_tables_to_csv(
    comp, tables_outputFlows, tables_inputFlows, MP_form_dict_reverse, size_dict, path
):
    # Extract input and output flow tables
    df1 = tables_outputFlows[comp]
    MP_size = []
    MP_form = []
    for x in df1.index:
        MP_size.append(size_dict[x[0]])
        MP_form.append(MP_form_dict_reverse[x[1:2]])
    df1.insert(0, "MP_size", MP_size)
    df1.insert(1, "MP_form", MP_form)

    df2 = tables_inputFlows[comp]
    MP_size = []
    MP_form = []
    for x in df2.index:
        MP_size.append(size_dict[x[0]])
        MP_form.append(MP_form_dict_reverse[x[1:2]])
    df2.insert(0, "MP_size", MP_size)
    df2.insert(1, "MP_form", MP_form)

    # Create a Pandas Excel writer using the file name and desired engine
    writer = pd.ExcelWriter(path + "/" + comp + "_mass_flows.xlsx", engine="xlsxwriter")

    # Write each dataframe to a different worksheet
    df1.to_excel(writer, sheet_name="Output_flows", index=False)
    df2.to_excel(writer, sheet_name="Input_flows", index=False)

    # Save the file
    writer.save()


def results_extended_by_compartment_to_csv(path, results_dict):
    # Create a Pandas Excel writer using the file name and desired engine
    writer = pd.ExcelWriter(path + "/" + "Results_SS.xlsx", engine="xlsxwriter")

    # Write each dataframe to a different worksheet (for each compartment)
    for comp in results_dict:
        df_results = results_dict[comp]
        df_results.to_excel(writer, sheet_name=comp, index=False)
    # Save the file
    writer.save()


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def plot_fractionDistribution_heatmap(Results_extended, fraction):
    # Define the order for the MP_Form labels
    mp_form_order = [
        "freeMP",
        "heterMP",
        "biofMP",
        "heterBiofMP",
    ]  # Replace with your desired order

    # Define the order for the Compartment labels
    compartment_order = [
        "Ocean_Surface_Water",
        "Ocean_Mixed_Water",
        "Ocean_Column_Water",
        "Coast_Surface_Water",
        "Coast_Column_Water",
        "Surface_Freshwater",
        "Bulk_Freshwater",
        "Sediment_Freshwater",
        "Sediment_Ocean",
        "Sediment_Coast",
        "Urban_Soil_Surface",
        "Urban_Soil",
        "Background_Soil_Surface",
        "Background_Soil",
        "Agricultural_Soil_Surface",
        "Agricultural_Soil",
        "Air",
    ]  # Replace with your desired order

    # Pivot the DataFrame to have one row per combination of MP_Form, Compartment, and Size_Fraction_um
    pivot_table = Results_extended.pivot_table(
        index=["MP_Form", "Size_Fraction_um"],
        columns="Compartment",
        values=fraction,
        aggfunc="mean",
    )

    # Reorder the rows based on mp_form_order and columns based on compartment_order
    pivot_table = pivot_table.loc[mp_form_order, compartment_order]

    # Apply log scale to the pivot table
    pivot_table_log = np.log10(pivot_table)

    # Replace -inf values with NaN
    pivot_table_log.replace(-np.inf, np.nan, inplace=True)

    # Stablish a lower limit
    # Set the lower limit for the values
    lower_limit = -14
    upper_limit = np.nanmax(pivot_table_log)

    # Replace values below the lower limit with NaN
    pivot_table_log = pivot_table_log.applymap(
        lambda x: np.nan if x < lower_limit else x
    )

    # Define a custom colormap with grey color for NaN values
    cmap = sns.color_palette("viridis", as_cmap=True)
    cmap.set_bad("white")

    # Plot the heatmap with logarithmic scale and custom colormap
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        pivot_table_log,
        cmap=cmap,
        cbar=True,
        cbar_kws={"label": "log10 (" + fraction + ") "},
        annot=False,
        linewidths=0.5,
        linecolor="grey",
        vmin=lower_limit,
        vmax=upper_limit,
    )

    # Set compartment labels to cover all size fractions underneath
    compartment_labels = pivot_table.columns
    compartment_label_positions = np.arange(len(compartment_labels)) + 0.5
    plt.xticks(
        ticks=compartment_label_positions, labels=compartment_labels, rotation=90
    )

    # Set MP_Form and Size_Fraction_um labels
    row_labels = [
        f"{mp_form} - {size_frac_um}" for mp_form, size_frac_um in pivot_table.index
    ]
    row_label_positions = np.arange(len(pivot_table.index)) + 0.5
    plt.yticks(ticks=row_label_positions, labels=row_labels, rotation=0)

    plt.title(
        "Heatmap of log10 ("
        + fraction
        + " by MP_Form, Compartment, and Size_Fraction_um"
    )
    plt.xlabel("Compartment", fontsize=14)
    plt.ylabel("MP_Form - Size_Fraction_um", fontsize=14)

    plt.tight_layout()
    plt.show()
