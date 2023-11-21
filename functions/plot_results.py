import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm

particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


def plot_bySize_total_number_particles(results_dict, comp_name, dict_size_coding):
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
    plt.close(fig)

    return fig


def plot_bySize_total_mass(results_dict, comp_name, dict_size_coding):
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

    plt.close(fig)

    return fig


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


def extract_output_table(comp, tables_outputFlows, MP_form_dict_reverse, size_dict):
    T = tables_outputFlows[comp]
    MP_size = []
    MP_form = []
    for x in T.index:
        MP_size.append(size_dict[x[0]])
        MP_form.append(MP_form_dict_reverse[x[1:2]])
    T.insert(0, "MP_size", MP_size)
    T.insert(1, "MP_form", MP_form)
    if sum(T.loc[:, T.columns[2:]].max()) == 0:
        print("All values are 0 and no heatmap can be printed")
    else:
        ax = plt.axes()
        sns.heatmap(
            T.loc[:, T.columns[2:]],
            xticklabels=True,
            yticklabels=True,
            norm=LogNorm(),
            linewidths=1,
            linecolor="grey",
            ax=ax,
        )
        plt.title("Output Flows for " + comp + " (g/s)")
        plt.xlabel("Process")
        plt.ylabel("Particle")
        plt.show()
    return T
