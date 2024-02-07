import numpy as np
import pandas as pd


def fillInteractions_fun_OOP(system_particle_object_list, SpeciesList, surfComp_list):
    # Asign loose rates
    elimination_rates = eliminationProcesses(system_particle_object_list, SpeciesList)

    interactions_df = pd.DataFrame(
        np.diag(elimination_rates), index=SpeciesList, columns=SpeciesList
    )

    # Asign interactions rates
    interactions_df_rows = []

    for sp1 in system_particle_object_list:
        interactions_df_rows.append(
            interactionProcess(
                sp1, interactions_df, system_particle_object_list, surfComp_list
            )
        )

    # interact3(sp1) for sp1 in interactions_df.index.to_list()]
    array = np.column_stack(
        interactions_df_rows
    )  # vstack it was set as column stack and was wrong!!
    interactions_df_sol = pd.DataFrame(array, index=SpeciesList, columns=SpeciesList)

    return interactions_df_sol.transpose()


def eliminationProcesses(system_particle_object_list, SpeciesList):
    # Estimate losses (diagonal):the diagonal of the dataframe corresponds to the losses of each species

    """create the array of values for the diagonal wich is the sum of all RC corresponding to one species:"""

    diag_list = []

    for sp in system_particle_object_list:
        dict_sp = sp.RateConstants

        # replace none values
        for k in dict_sp.keys():
            if dict_sp[k] == None:
                dict_sp[k] = 0
            else:
                pass

        losses = []
        for (
            i
        ) in (
            dict_sp.keys()
        ):  # Do we need to remove fragemtation of smallest size bin and only take into consideration desintegration??
            if i == "k_fragmentation":
                if type(dict_sp[i]) == tuple:
                    losses.append(dict_sp[i][0])
                else:
                    losses.append(dict_sp[i])

            elif type(dict_sp[i]) == tuple or type(dict_sp[i]) == list:
                losses.append(sum(dict_sp[i]))
            else:
                losses.append(dict_sp[i])

        diag_list.append(-(sum(losses)))

    return diag_list


def inboxProcess(sp1, sp2, surfComp_list):
    # If same compartment (compartment processes)
    if sp1.Pcode[2:] == sp2.Pcode[2:]:
        # Only different size bins --> Fragmentation

        if sp1.Pcode[1:] == sp2.Pcode[1:] and sp1.Pcode[0] != sp2.Pcode[0]:

            # We reformulate fractionation as it is not a happening only for consecutive size Bins(bigger to next smaller) but using the fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html)
            #  fsd = is a 2D matrix, but as fragmentation only happens from larger to smaller size classes, the matrix is zero everywhere except its lower left triangle.

            # Set fsd. This is an example of distribution but can be changed and also subject to beta factor as described in the link avobe. To be discussed
            fsd = np.array(
                [
                    [0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0],
                    [0.5, 0.5, 0, 0, 0],
                    [0.6, 0.2, 0.2, 0, 0],
                    [0.7, 0.15, 0.1, 0.05, 0],
                ]
            )
            # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class

            size_dict = {chr(i): i - ord("a") for i in range(ord("a"), ord("e") + 1)}

            fsd_index_i = size_dict[sp2.Pcode[0]]
            fsd_index_j = size_dict[sp1.Pcode[0]]
            frag_factor = fsd[fsd_index_i, fsd_index_j]
            if type(sp2.RateConstants["k_fragmentation"]) is tuple:
                frag = sp2.RateConstants["k_fragmentation"]

                sol = frag[0] * frag_factor
            else:
                sol = sp2.RateConstants["k_fragmentation"] * frag_factor

        # Different aggergation states (same size)--> heteroagg, biofouling,defouling and agg-breackup

        elif sp1.Pcode[0] == sp2.Pcode[0] and sp1.Pcode[1] != sp2.Pcode[1]:
            # heteroaggregation from A-->B or from C-->D
            if (sp2.Pcode[1] == "A" and sp1.Pcode[1] == "B") or (
                sp2.Pcode[1] == "C" and sp1.Pcode[1] == "D"
            ):
                process = "heteroaggregation"
                if process in sp2.Pcompartment.processess:
                    sol = sp2.RateConstants["k_" + process]
                else:
                    sol = 0

            # heteroaggregate breackup from B-->A and from D-->C
            elif (sp2.Pcode[1] == "B" and sp1.Pcode[1] == "A") or (
                sp2.Pcode[1] == "D" and sp1.Pcode[1] == "C"
            ):
                process = "heteroaggregate_breackup"
                if process in sp2.Pcompartment.processess:
                    sol = sp2.RateConstants["k_" + process]
                else:
                    sol = 0

            # Biofouling from A-->C or from B-->D
            elif (sp2.Pcode[1] == "A" and sp1.Pcode[1] == "C") or (
                sp2.Pcode[1] == "B" and sp1.Pcode[1] == "D"
            ):
                process = "biofouling"
                if process in sp2.Pcompartment.processess:
                    sol = sp2.RateConstants["k_" + process]
                else:
                    sol = 0

            # Defouling from C-->A or from D-->B
            elif (sp2.Pcode[1] == "C" and sp1.Pcode[1] == "A") or (
                sp2.Pcode[1] == "D" and sp1.Pcode[1] == "B"
            ):
                process = "defouling"
                if process in sp2.Pcompartment.processess:
                    sol = sp2.RateConstants["k_" + process]
                else:
                    sol = 0

            else:
                sol = 0
        else:
            sol = 0

    # Different compartments--> Transport processess
    # settling, rising, mixing, resusp, advective transport, difussion, runoff, percolation?

    # if compartments are in the list of compartment connexions
    # check if same agg form and size to select process of connexion for compartment and assign rate constant, else process has rate of 0

    elif sp1.Pcompartment.Cname in sp2.Pcompartment.connexions:
        # transport between compartments only for same aggregation state and same particle size
        if sp1.Pcode[:2] == sp2.Pcode[:2]:
            process = sp2.Pcompartment.connexions[sp1.Pcompartment.Cname]
            if type(process) == list:
                sol2 = []
                for p in process:
                    if p == "dry_depossition":
                        element = sp1.Pcompartment.Cname
                        # Iterate over the list using enumerate()
                        for index, value in enumerate(surfComp_list):
                            if value == element:
                                position = index
                        sol2.append(sp2.RateConstants["k_" + p][position])

                    else:
                        sol2.append(sp2.RateConstants["k_" + p])
                sol = sum(sol2)
            else:
                if process == "runoff_transport":
                    # Runoff can happen from soil surface to freshwater surface and to coast surface water and it will happen in different proportions. we stablish the fraction that goes into each water body through fdd : the matrix containing the fractions of the runoff of each surface soil type that goes into each surface water body (only coast and freshwater)
                    ro_comp = {
                        "Urban_Soil_Surface": 0,
                        "Background_Soil_Surface": 1,
                        "Agricultural_Soil_Surface": 2,
                    }
                    recieving_comp = {"Coast_Surface_Water": 0, "Surface_Freshwater": 1}

                    # In this example of fdd all runoff goes to surface freshwater. To be discussed later
                    fro = np.array([[0, 1], [0, 1], [0, 1]])

                    sol = (
                        sp2.RateConstants["k_" + process]
                        * fro[
                            ro_comp[sp2.Pcompartment.Cname],
                            recieving_comp[sp1.Pcompartment.Cname],
                        ]
                    )
                else:
                    sol = sp2.RateConstants["k_" + process]
        else:
            sol = 0
    else:
        sol = 0

    return sol


def interactionProcess(
    sp1, interactions_df, system_particle_object_list, surfComp_list
):
    sol = []
    for sp2 in system_particle_object_list:
        # Same particle in the same box and compartment (losses)
        if sp1.Pcode == sp2.Pcode:
            sol.append(interactions_df[sp2.Pcode][sp1.Pcode])

        # Different particle or different river section or compartment
        else:
            # Same box (i.e. river section RS)--> In box processes

            if sp1.Pcode.split("_")[1] == sp2.Pcode.split("_")[1]:
                sol.append(inboxProcess(sp1, sp2, surfComp_list))

            # Different Box but same particle in same compartment (Full Multi version where more than 1 box (i.e. river sections)) -->Transport (advection or sediment transport determined by flow_connectivity file)

            elif sp2.Pcode.split("_")[0] == sp1.Pcode.split("_")[0]:
                sol.append(transportProcess(sp1, sp2, river_flows))
                # """Pending work on transport for The Full Multi"""
            else:
                sol.append(0)

    return sol


def transportProcess(sp1, sp2, RC_df, river_flows):
    J = int(sp1[:-3]) + 1
    I = int(sp2[:-3]) + 1
    flowI_df = river_flows[river_flows.Region_I == I]
    if J in flowI_df.Region_J.tolist():
        if sp1[-3] != "4":
            if isinstance(RC_df[sp2]["advection"], (int, float)):
                solution = RC_df[sp2]["advection"]
            else:
                idx_ad = np.where(flowI_df.Region_J == J)[0][0]
                solution = RC_df[sp2]["advection"][idx_ad]
        else:
            solution = RC_df[sp2]["sedTransport"]
    else:
        solution = 0

    return solution
