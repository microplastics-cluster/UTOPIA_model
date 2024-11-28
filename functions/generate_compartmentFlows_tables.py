import pandas as pd


def estimate_outFlows(system_particle_object_list, dict_comp):
    # Estimate mass outflows
    for p in system_particle_object_list:
        p.outFlow_mass_g_s = {}
        p.outFlow_number_g_s = {}
        for c in p.RateConstants:
            if type(p.RateConstants[c]) == list:
                p.outFlow_mass_g_s[c] = [R * p.Pmass_g_SS for R in p.RateConstants[c]]
                p.outFlow_number_g_s[c] = [R * p.Pnum_SS for R in p.RateConstants[c]]
            else:
                p.outFlow_mass_g_s[c] = p.RateConstants[c] * p.Pmass_g_SS
                p.outFlow_number_g_s[c] = p.RateConstants[c] * p.Pnum_SS

    # Tables of output flows per compartmet
    tables_outputFlows_mass = {}
    tables_outputFlows_number = {}
    for c in list(dict_comp.keys()):
        part_dic_mass = {}
        part_dic_number = {}
        for p in system_particle_object_list:
            if p.Pcompartment.Cname == c:
                part_dic_mass[p.Pcode] = pd.DataFrame.from_dict(
                    p.outFlow_mass_g_s, orient="index"
                )
                part_dic_number[p.Pcode] = pd.DataFrame.from_dict(
                    p.outFlow_number_g_s, orient="index"
                )
        tables_outputFlows_mass[c] = pd.concat(part_dic_mass, axis=1).transpose()
        tables_outputFlows_number[c] = pd.concat(part_dic_number, axis=1).transpose()

    for k in tables_outputFlows_mass:
        tables_outputFlows_mass[k] = (
            tables_outputFlows_mass[k].reset_index(level=1).drop("level_1", axis=1)
        )
        tables_outputFlows_number[k] = (
            tables_outputFlows_number[k].reset_index(level=1).drop("level_1", axis=1)
        )

    return tables_outputFlows_mass, tables_outputFlows_number


# Estimate imput flows from transport from other compartments


def estimate_inFlows(
    tables_outputFlows, tables_outputFlows_number, dict_comp, surfComp_list
):
    ##Tables of recieving flows through transport from other compartments
    tables_inputFlows = {}
    tables_inputFlows_num = {}
    for comp in list(dict_comp.keys()):
        comp_input_flows = []
        comp_input_flows_num = []
        for e_comp in dict_comp:
            if comp in dict_comp[e_comp].connexions:
                inpProc = dict_comp[e_comp].connexions[comp]
                if (
                    type(inpProc) == list
                ):  # When there is more than one process of inflow into the compartment
                    df_inflows = tables_outputFlows[e_comp].loc[
                        :, ["k_" + ele for ele in inpProc]
                    ]
                    df_inflows_num = tables_outputFlows_number[e_comp].loc[
                        :, ["k_" + ele for ele in inpProc]
                    ]

                    for proc in inpProc:
                        if proc == "dry_deposition":
                            position = surfComp_list.index(comp)
                            df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                lambda x: x[position] if isinstance(x, list) else x
                            )
                            df_inflows_num["k_" + proc] = df_inflows_num[
                                "k_" + proc
                            ].apply(lambda x: x[position] if isinstance(x, list) else x)

                        elif proc == "mixing":

                            if (
                                e_comp == "Ocean_Mixed_Water"
                                and comp == "Ocean_Surface_Water"
                            ):
                                df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                    lambda x: x[0] if isinstance(x, list) else x
                                )
                                df_inflows_num["k_" + proc] = df_inflows_num[
                                    "k_" + proc
                                ].apply(lambda x: x[0] if isinstance(x, list) else x)

                            elif (
                                e_comp == "Ocean_Mixed_Water"
                                and comp == "Ocean_Column_Water"
                            ):
                                df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                    lambda x: x[1] if isinstance(x, list) else x
                                )
                                df_inflows_num["k_" + proc] = df_inflows_num[
                                    "k_" + proc
                                ].apply(lambda x: x[1] if isinstance(x, list) else x)
                            else:
                                pass
                            # Revisit for percollation and tillage
                        else:
                            pass
                    comp_input_flows.append(df_inflows)
                    comp_input_flows_num.append(df_inflows_num)

                else:
                    df_inflows = (
                        tables_outputFlows[e_comp].loc[:, "k_" + inpProc].to_frame()
                    )
                    df_inflows_num = (
                        tables_outputFlows_number[e_comp]
                        .loc[:, "k_" + inpProc]
                        .to_frame()
                    )
                    for ele in df_inflows["k_" + inpProc]:
                        if type(ele) == list:
                            connecting_comp = {
                                key: value
                                for key, value in dict_comp[e_comp].connexions.items()
                                if value == inpProc
                            }
                            poss_dict = {
                                key: index
                                for index, key in enumerate(connecting_comp.keys())
                            }
                            possition = poss_dict[comp]
                            df_inflows["k_" + inpProc] = df_inflows[
                                "k_" + inpProc
                            ].apply(
                                lambda x: x[possition] if isinstance(x, list) else x
                            )
                            df_inflows_num["k_" + inpProc] = df_inflows_num[
                                "k_" + inpProc
                            ].apply(
                                lambda x: x[possition] if isinstance(x, list) else x
                            )

                        else:
                            pass
                    comp_input_flows.append(df_inflows)
                    comp_input_flows_num.append(df_inflows_num)
            else:
                pass

        tables_inputFlows[comp] = pd.concat(comp_input_flows).fillna(0)
        tables_inputFlows_num[comp] = pd.concat(comp_input_flows_num).fillna(0)
    return tables_inputFlows, tables_inputFlows_num


def generate_flows_dict(
    tables_outputFlows, tables_inputFlows, size_dict, MP_form_dict_reverse
):
    flows_dict = dict()
    flows_dict["input_flows"] = {}
    flows_dict["output_flows"] = {}

    # Decode index in input and output flow tables
    for comp in tables_outputFlows.keys():
        df1 = tables_outputFlows[comp].copy()
        MP_size_df1 = []
        MP_form_df1 = []
        for x in df1.index:
            MP_size_df1.append(size_dict[x[0]])
            MP_form_df1.append(MP_form_dict_reverse[x[1:2]])

        df1.insert(0, "MP_size", MP_size_df1)
        df1.insert(1, "MP_form", MP_form_df1)
        flows_dict["output_flows"][comp] = df1

    for comp in tables_inputFlows:
        df2 = tables_inputFlows[comp].copy()
        MP_size_df2 = []
        MP_form_df2 = []
        for y in df2.index:
            MP_size_df2.append(size_dict[y[0]])
            MP_form_df2.append(MP_form_dict_reverse[y[1:2]])
        df2.insert(0, "MP_size", MP_size_df2)
        df2.insert(1, "MP_form", MP_form_df2)
        flows_dict["input_flows"][comp] = df2

    return flows_dict


# Function to handle summing lists and individual elements
def handle_value(value):
    if isinstance(value, list):
        return sum(value)
    return value


def sum_column_values(column):
    return sum(handle_value(value) for value in column)


def process_flows(compartment, size_fraction, mp_form, flow_type, flows_dict):
    """Process flows (inflows or outflows) for a given compartment, size fraction, and MP form."""
    df_comp = flows_dict[flow_type][compartment]
    df_filtered = df_comp[
        (df_comp["MP_form"] == mp_form) & (df_comp["MP_size"] == size_fraction)
    ]
    df_cleaned = df_filtered.drop(["MP_size", "MP_form"], axis=1)
    return {col: sum_column_values(df_cleaned[col]) for col in df_cleaned.columns}


def process_flows_comp(compartment, flow_type, flows_dict):
    """Process flows (inflows or outflows) for a given compartment"""
    df_comp = flows_dict[flow_type][compartment]
    df_cleaned = df_comp.drop(["MP_size", "MP_form"], axis=1)

    return {col: sum_column_values(df_cleaned[col]) for col in df_cleaned.columns}


def addFlows_to_results_df(Results_extended, flows_dict_mass, flows_dict_num):
    """Calculate inflows and outflows (mass and number) and update Results_extended."""
    inflows_mass_list = []
    inflows_num_list = []
    outflows_mass_list = []
    outflows_num_list = []

    for n in range(len(Results_extended)):
        compartment = Results_extended.iloc[n]["Compartment"]
        size_fraction = Results_extended.iloc[n]["Size_Fraction_um"]
        mp_form = Results_extended.iloc[n]["MP_Form"]

        # Calculate inflows and outflows for mass
        inflows_mass = process_flows(
            compartment, size_fraction, mp_form, "input_flows", flows_dict_mass
        )
        outflows_mass = process_flows(
            compartment, size_fraction, mp_form, "output_flows", flows_dict_mass
        )
        inflows_mass_list.append(inflows_mass)
        outflows_mass_list.append(outflows_mass)

        # Calculate inflows and outflows for number
        inflows_num = process_flows(
            compartment, size_fraction, mp_form, "input_flows", flows_dict_num
        )
        outflows_num = process_flows(
            compartment, size_fraction, mp_form, "output_flows", flows_dict_num
        )
        inflows_num_list.append(inflows_num)
        outflows_num_list.append(outflows_num)

    # Update the Results_extended DataFrame with the calculated flows
    Results_extended["inflows_g_s"] = inflows_mass_list
    Results_extended["inflows_num_s"] = inflows_num_list
    Results_extended["outflows_g_s"] = outflows_mass_list
    Results_extended["outflows_num_s"] = outflows_num_list

    return Results_extended


def addFlows_to_results_df_comp(mass_dist_comp, flows_dict_mass, flows_dict_num):
    """Calculate inflows and outflows (mass and number) and update Results_extended."""
    inflows_mass_list = []
    inflows_num_list = []
    outflows_mass_list = []
    outflows_num_list = []

    for n in range(len(mass_dist_comp)):
        compartment = mass_dist_comp.iloc[n]["Compartments"]

        # Calculate inflows and outflows for mass
        inflows_mass = process_flows_comp(compartment, "input_flows", flows_dict_mass)
        outflows_mass = process_flows_comp(compartment, "output_flows", flows_dict_mass)
        inflows_mass_list.append(inflows_mass)
        outflows_mass_list.append(outflows_mass)

        # Calculate inflows and outflows for number
        inflows_num = process_flows_comp(compartment, "input_flows", flows_dict_num)
        outflows_num = process_flows_comp(compartment, "output_flows", flows_dict_num)
        inflows_num_list.append(inflows_num)
        outflows_num_list.append(outflows_num)

    # Update the Results_extended DataFrame with the calculated flows
    mass_dist_comp["inflows_g_s"] = inflows_mass_list
    mass_dist_comp["inflows_num_s"] = inflows_num_list
    mass_dist_comp["outflows_g_s"] = outflows_mass_list
    mass_dist_comp["outflows_num_s"] = outflows_num_list

    return mass_dist_comp


def add_output_flow_conexions(
    results_by_comp,
    dict_comp,
    outputflow_type="outflows_g_s",
    inputflow_type="inflows_g_s",
):
    # Genrate table of output flows connexions beteen compartments to b used in the mass flow diagram (ti me developed)

    outflow_conexions_g_s = []
    for c in results_by_comp["Compartments"]:
        conexions = dict_comp[c].connexions
        outfows = results_by_comp[results_by_comp["Compartments"] == c][
            outputflow_type
        ].values[0]
        outflow_conexions = {}
        for key, value in conexions.items():
            # Check if the value is a list or a single string
            if isinstance(value, list):
                # Create a dictionary for each element in the list
                if c == "Ocean_Mixed_Water":
                    inflows_col = results_by_comp[
                        results_by_comp["Compartments"] == key
                    ][inputflow_type].values[0]
                    outflow_conexions[key] = {
                        item: (
                            inflows_col["k_" + item]
                            if item == "mixing"
                            else outfows["k_" + item]
                        )
                        for item in value
                    }
                else:
                    outflow_conexions[key] = {
                        item: outfows["k_" + item] for item in value
                    }
            else:
                # If it's a single string, directly map it
                outflow_conexions[key] = {value: outfows["k_" + value]}

        # Print the modified dictionary
        outflow_conexions_g_s.append(outflow_conexions)
    results_by_comp["outflow_conexions_" + outputflow_type[9:]] = outflow_conexions_g_s
    return results_by_comp
