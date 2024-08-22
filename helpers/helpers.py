import matplotlib.pyplot as plt


def generate_system_species_list(
    system_particle_object_list, MPforms_list, compartmentNames_list, boxNames_list
):
    particle_sizes_coding = {"mp1": "a", "mp2": "b", "mp3": "c", "mp4": "d", "mp5": "e"}

    particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))

    particle_compartmentCoding = dict(
        zip(compartmentNames_list, list(range(len(compartmentNames_list))))
    )

    def particle_nameCoding(particle, boxNames_list):
        # if len(boxNames_list) != 1:

        particle_sizeCode = particle_sizes_coding[particle.Pname[0:3]]
        particle_formCode = particle_forms_coding[particle.Pform]
        particle_compartmentCode = particle_compartmentCoding[
            particle.Pcompartment.Cname
        ]
        particle_boxCode = particle.Pcompartment.CBox.Bname

        particleCode = (
            particle_sizeCode
            + particle_formCode
            + str(particle_compartmentCode)
            + "_"
            + particle_boxCode
        )
        # else:
        #     particle_sizeCode = particle_sizes_coding[particle.Pname[0:3]]
        #     particle_formCode = particle_forms_coding[particle.Pform]
        #     particle_compartmentCode = particle_compartmentCoding[
        #         particle.Pcompartment.Cname
        #     ]

        #     particleCode = (
        #         particle_sizeCode + particle_formCode + str(particle_compartmentCode)
        #     )

        return particleCode

    SpeciesList = []
    for particle in system_particle_object_list:
        SpeciesList.append(particle_nameCoding(particle, boxNames_list))
        particle.Pcode = particle_nameCoding(particle, boxNames_list)

    return SpeciesList


# Plot rate constat values for comparison
def plot_rate_constants(RC_df):
    processList = processList = [k for k in RC_df.columns if "k" in k]

    df_RC = RC_df[processList]

    fig = df_RC.plot(
        title="Rate constant values (s-1)",
        subplots=True,
        figsize=(10, 15),
        sharex=True,
        fontsize=12,
        stacked=True,
    )
    plt.savefig("rateConstants.png")

    plt.show()


import numpy as np
import pandas as pd


def timeLimit_particles_RC(system_particle_object_list, lim):
    for particle in system_particle_object_list:
        for k in particle.RateConstants:
            if (
                particle.RateConstants[k] is not None
                and particle.RateConstants[k] > lim
            ):
                particle.RateConstants[k] = lim
            else:
                pass
    return system_particle_object_list


def change_RC_value(system_particle_object_list, rc_name, rc_val):
    for p in system_particle_object_list:
        for k in p.RateConstants:
            if k == rc_name:
                p.RateConstants[k] = rc_val
            else:
                pass
    return system_particle_object_list


# function to convert mass to number
def mass_to_num(mass_g, volume_m3, density_kg_m3):
    number = mass_g / 1000 / density_kg_m3 / volume_m3
    # number of particles has always to be integer?
    return number


# function to convert number to mass
def num_to_mass(number, volume_m3, density_kg_m3):
    mass_g = number * volume_m3 * density_kg_m3 * 1000
    return mass_g


def extract_inflows_outflows(flows_dict_mass, comp, MP_form, MP_size):

    # function to extract the input and output flows of the system (cell selection)

    # extract inflows
    df_i = flows_dict_mass["input_flows"][comp]
    df_ii = df_i.loc[(df_i["MP_form"] == MP_form) & (df_i["MP_size"] == MP_size)]
    df_iii = df_ii.drop(["MP_size", "MP_form"], axis=1)
    list_iflows = [k for k in df_iii]
    list_inflow_val = [sum(df_iii[col]) for col in list_iflows]
    inflow_p = [round((v / sum(list_inflow_val)) * 100, 4) for v in list_inflow_val]

    pd_inputFlows = pd.DataFrame(
        {"Inflows": list_iflows, "Rate_g_s": list_inflow_val, "%": inflow_p}
    )

    # extract outflows
    df_o = flows_dict_mass["output_flows"][comp]
    df_oo = df_o.loc[(df_o["MP_form"] == MP_form) & (df_o["MP_size"] == MP_size)]
    df_oo.reset_index(drop=True, inplace=True)
    df_ooo = df_oo.drop(["MP_size", "MP_form"], axis=1)

    outflow_val = [
        sum(df_ooo[ko]) if type(df_ooo[ko][0]) != list else sum(df_ooo[ko][0])
        for ko in df_ooo
    ]
    df_ooo.loc[0] = outflow_val
    list_outflows = [ko for ko in df_ooo.columns if (df_ooo[ko] != 0).any()]
    list_outflow_val = [
        value for column in df_ooo.columns for value in df_ooo[column] if value != 0
    ]
    outflow_p = [round((v / sum(list_outflow_val)) * 100, 4) for v in list_outflow_val]

    pd_outflows = dict(zip(list_outflows, (list_outflow_val, outflow_p)))
    # pd.DataFrame(
    #     {"Outflows": list_outflows, "Rate_g_s": list_outflow_val, "%": outflow_p}
    # )

    return pd_inputFlows, pd_outflows


def extract_inflows_outflows_comp(flows_dict_mass, comp):

    # function to extract the input and output flows of the system (per compartmet)

    # extract inflows
    df_i = flows_dict_mass["input_flows"][comp]
    df_ii = df_i.drop(["MP_size", "MP_form"], axis=1)
    list_iflows = [k for k in df_ii]
    list_iflows
    list_inflow_val = [sum(df_ii[col]) for col in list_iflows]
    inflow_p = [round((v / sum(list_inflow_val)) * 100, 4) for v in list_inflow_val]
    pd_inputFlows = pd.DataFrame(
        {"Inflows": list_iflows, "Rate_g_s": list_inflow_val, "%": inflow_p}
    )
    pd_inputFlows = pd.DataFrame(
        {"Inflows": list_iflows, "Rate_g_s": list_inflow_val, "%": inflow_p}
    )

    # extract outflows
    df_o = flows_dict_mass["output_flows"][comp]
    df_o.reset_index(drop=True, inplace=True)
    df_oo = df_o.drop(["MP_size", "MP_form"], axis=1)

    outflow_val = [
        sum(df_oo[ko]) if type(df_oo[ko][0]) != list else sum(df_oo[ko][0])
        for ko in df_oo
    ]
    list_outflows = [
        df_oo.columns[i] for i in range(len(df_oo.columns)) if outflow_val[i] != 0
    ]
    list_outflow_val = [val for val in outflow_val if val != 0]

    outflow_p = [round((v / sum(list_outflow_val)) * 100, 4) for v in list_outflow_val]

    pd_outflows = pd.DataFrame(
        {"Outflows": list_outflows, "Rate_g_s": list_outflow_val, "%": outflow_p}
    )
    pd_outflows

    return pd_inputFlows, pd_outflows
