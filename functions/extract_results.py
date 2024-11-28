# function to organise results into a dictionary of compartments, each compartment containing a dictionary of aggregation states, of wich contain number of particles per size fraction
import pandas as pd


def extract_by_comp(Results, particle_compartmentCoding):
    Results_comp_dict = {}
    for comp in particle_compartmentCoding:
        key = comp
        values = []
        for s in Results["species"]:
            if str(particle_compartmentCoding[comp]) == s[2:-7]:
                values.append(s)

            else:
                pass
        Results_comp_dict[key] = Results[Results["species"].isin(values)]

    return Results_comp_dict


def extract_by_aggSt(Results_comp_dict, particle_forms_coding):
    Results_comp_organiced = {}
    for comp in Results_comp_dict:
        Results_aggSt_dict = {}
        key1 = comp
        df = Results_comp_dict[comp]
        for aggst in particle_forms_coding:
            key2 = aggst
            values = []
            for c in df["species"]:
                if particle_forms_coding[aggst] == c[1]:
                    values.append(c)
                else:
                    pass
            Results_aggSt_dict[key2] = df[df["species"].isin(values)]
        Results_comp_organiced[key1] = Results_aggSt_dict
    return Results_comp_organiced


def extract_by_size(Results):
    sizes = ["a", "b", "c", "d", "e"]
    Results_size_dict = {}
    for x in sizes:
        key = x
        values = []
        for s in Results["species"]:
            if x in s[0:2]:
                values.append(s)
        Results_size_dict[key] = Results[Results["species"].isin(values)]
    print(Results_size_dict)
    return Results_size_dict


def estimate_fractions(Results):
    total_mass = sum(Results["mass_g"])
    total_number = sum(Results["number_of_particles"])
    Results_extended = Results.copy()
    Results_extended.loc[:, "mass_fraction"] = [
        x / total_mass for x in Results["mass_g"]
    ]
    Results_extended.loc[:, "number_fraction"] = [
        x / total_number for x in Results["number_of_particles"]
    ]

    mass_fraction_df = Results_extended.loc[
        :, ["Compartment", "MP_Form", "Size_Fraction_um", "mass_fraction"]
    ]

    number_fraction_df = Results_extended.loc[
        :, ["Compartment", "MP_Form", "Size_Fraction_um", "number_fraction"]
    ]

    # Short values (descending)
    mf_shorted = mass_fraction_df.sort_values("mass_fraction", ascending=False)

    nf_shorted = number_fraction_df.sort_values("number_fraction", ascending=False)

    # print(mf_shorted[:10])
    # print(nf_shorted[:10])

    return Results_extended, mf_shorted, nf_shorted


def results_subset(Results_df, crit, cond):
    # if type(cond) == list:
    #     out_df={}
    #     for ele in cond:
    #         out_df[ele]=Results_df[Results_df[crit]==ele]
    #     else:
    out_df = Results_df[Results_df[crit] == cond]
    return out_df


def extract_results_by_compartment(Results_extended, dict_comp):
    mass_g = []
    particle_number = []
    mass_frac_100 = []
    num_frac_100 = []
    mass_conc_g_m3 = []
    num_conc = []
    for comp in list(dict_comp.keys()):
        mass_g.append(
            sum(Results_extended[Results_extended["Compartment"] == comp]["mass_g"])
        )
        particle_number.append(
            sum(
                Results_extended[Results_extended["Compartment"] == comp][
                    "number_of_particles"
                ]
            )
        )
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
        mass_conc_g_m3.append(
            sum(
                Results_extended[Results_extended["Compartment"] == comp][
                    "concentration_g_m3"
                ]
            )
        )
        num_conc.append(
            sum(
                Results_extended[Results_extended["Compartment"] == comp][
                    "concentration_num_m3"
                ]
            )
        )

    mass_dist_comp = pd.DataFrame(columns=["Compartments"])
    mass_dist_comp["Compartments"] = list(dict_comp.keys())
    mass_dist_comp["mass_g"] = mass_g
    mass_dist_comp["number_of_particles"] = particle_number
    mass_dist_comp["%_mass"] = mass_frac_100
    mass_dist_comp["%_number"] = num_frac_100
    mass_dist_comp["Concentration_g_m3"] = mass_conc_g_m3
    mass_dist_comp["Concentration_num_m3"] = num_conc

    return mass_dist_comp
