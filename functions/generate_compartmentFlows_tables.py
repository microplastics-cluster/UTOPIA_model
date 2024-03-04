import pandas as pd


def estimate_outFlows(system_particle_object_list, dict_comp):
    # Estimate outflows
    for p in system_particle_object_list:
        p.outFlow_mass_g_s = {}
        for c in p.RateConstants:
            if type(p.RateConstants[c]) == list:
                p.outFlow_mass_g_s[c] = [R * p.Pmass_g_SS for R in p.RateConstants[c]]
            else:
                p.outFlow_mass_g_s[c] = p.RateConstants[c] * p.Pmass_g_SS

    # Tables of output flows per compartmet
    tables_outputFlows = {}
    for c in list(dict_comp.keys()):
        part_dic = {}
        for p in system_particle_object_list:
            if p.Pcompartment.Cname == c:
                part_dic[p.Pcode] = pd.DataFrame.from_dict(
                    p.outFlow_mass_g_s, orient="index"
                )
        tables_outputFlows[c] = pd.concat(part_dic, axis=1).transpose()
    for k in tables_outputFlows:
        tables_outputFlows[k] = (
            tables_outputFlows[k].reset_index(level=1).drop("level_1", axis=1)
        )
    return tables_outputFlows


# Estimate imput flows from transport from other compartments


def estimate_inFlows(tables_outputFlows, dict_comp, surfComp_list):
    ##Tables of recieving flows through transport from other compartments
    tables_inputFlows = {}
    for comp in list(dict_comp.keys()):
        comp_input_flows = []
        for e_comp in dict_comp:
            if comp in dict_comp[e_comp].connexions:
                inpProc = dict_comp[e_comp].connexions[comp]
                if (
                    type(inpProc) == list
                ):  # When there is more than one process of inflow into the compartment
                    df_inflows = tables_outputFlows[e_comp].loc[
                        :, ["k_" + ele for ele in inpProc]
                    ]

                    for proc in inpProc:
                        if proc == "dry_depossition":
                            position = surfComp_list.index(comp)
                            df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                lambda x: x[position] if isinstance(x, list) else x
                            )
                        elif proc == "mixing":

                            if (
                                e_comp == "Ocean_Mixed_Water"
                                and comp == "Ocean_Surface_Water"
                            ):
                                df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                    lambda x: x[0] if isinstance(x, list) else x
                                )
                            elif (
                                e_comp == "Ocean_Mixed_Water"
                                and comp == "Ocean_Column_Water"
                            ):
                                df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                    lambda x: x[1] if isinstance(x, list) else x
                                )
                            else:
                                pass
                            # Revisit for percollation and tillage
                        else:
                            pass
                    comp_input_flows.append(df_inflows)
                else:
                    df_inflows = (
                        tables_outputFlows[e_comp].loc[:, "k_" + inpProc].to_frame()
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

                        else:
                            pass
                    comp_input_flows.append(df_inflows)
            else:
                pass

        tables_inputFlows[comp] = pd.concat(comp_input_flows).fillna(0)
    return tables_inputFlows