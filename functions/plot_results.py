import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
import os

import pandas as pd

particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


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


def plot_heatmap_RC(comp, rateConstants_df, path_RC):
    T = rateConstants_df[rateConstants_df["Compartment"] == comp]

    selected_columns = T.columns[4:]
    column1 = T["MP_form"]
    column2 = T["Size_Bin"]

    # Select the desired columns
    df_selected = T[selected_columns]

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
    writer = pd.ExcelWriter(comp + "_mass_flows.xlsx", engine="xlsxwriter")

    # Write each dataframe to a different worksheet
    df1.to_excel(writer, sheet_name="Output_flows", index=False)
    df2.to_excel(writer, sheet_name="Input_flows", index=False)

    # Save the file
    writer.save()
