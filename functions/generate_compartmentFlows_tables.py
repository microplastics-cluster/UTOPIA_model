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
                        if proc == "dry_depossition":
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
