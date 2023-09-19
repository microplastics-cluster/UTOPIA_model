import numpy as np
import pandas as pd

def fillInteractions_Knames(system_particle_object_list, SpeciesList):

    # Asign loose rates
    elimination_rates = eliminationProcesses_knames(system_particle_object_list, SpeciesList)

    interactions_df_knames = pd.DataFrame(
        np.diag(elimination_rates), index=SpeciesList, columns=SpeciesList
    )

    # Asign interactions rates
    interactions_df_rows = []

    for sp1 in system_particle_object_list:
        interactions_df_rows.append(
            interactionProcess_knames(sp1, interactions_df_knames, system_particle_object_list)
        )

    # interact3(sp1) for sp1 in interactions_df.index.to_list()]
    array = np.column_stack(
        interactions_df_rows
    )  # vstack it was set as column stack and was wrong!!
    interactions_df_sol = pd.DataFrame(array, index=SpeciesList, columns=SpeciesList)

    return interactions_df_sol


def eliminationProcesses_knames(system_particle_object_list, SpeciesList):
    # Estimate losses (diagonal):the diagonal of the dataframe corresponds to the losses of each species

    """create the array of values for the diagonal wich is the sum of all RC corresponding to one species:"""

    diag_list = []
    for sp in system_particle_object_list:
        dict_sp=sp.RateConstants
        losses = []
        for i in dict_sp.keys():
            if dict_sp[i]!= 0:
                losses.append(i)
            else:
                pass
        

        diag_list.append("+".join(losses))

    """Revisit losses!!"""
    return diag_list


def inboxProcess_knames(sp1, sp2):

    # If same compartment (compartment processes)
    if sp1.Pcode[2:] == sp2.Pcode[2:]:

        # Only different size bins --> Fragmentation

        if sp1.Pcode[1:] == sp2.Pcode[1:] and sp1.Pcode[0] != sp2.Pcode[0]:

            # fragmentation only will occur from bigger to smaller and in consecutive sieBins Sizebin = sp[-3]
            if (
                (sp2.Pcode[0] == "b" and sp1.Pcode[0] == "a")
                or (sp2.Pcode[0] == "c" and sp1.Pcode[0] == "b")
                or (sp2.Pcode[0] == "d" and sp1.Pcode[0] == "c")
                or (sp2.Pcode[0] == "e" and sp1.Pcode[0] == "d")
            ):

                sol = "k_fragmentation"
            else:
                sol=0

        # Different aggergation states (same size)--> heteroagg, biofouling,defouling and agg-breackup

        elif sp1.Pcode[0] == sp2.Pcode[0] and sp1.Pcode[1] != sp2.Pcode[1]:

            # heteroaggregation from A-->B or from C-->D
            if (sp2.Pcode[1] == "A" and sp1.Pcode[1] == "B") or (
                sp2.Pcode[1] == "C" and sp1.Pcode[1] == "D"
            ):
                process = "heteroaggregation"
                if process in sp2.Pcompartment.processess:
                    sol = "k_" + process
                else:
                    sol=0

            # heteroaggregate breackup from B-->A and from D-->C
            elif (sp2.Pcode[1] == "B" and sp1.Pcode[1] == "A") or (
                sp2.Pcode[1] == "D" and sp1.Pcode[1] == "C"
            ):
                process = "heteroaggregate_breackup"
                if process in sp2.Pcompartment.processess:
                    sol = "k_" + process
                else:
                    sol=0

            # Biofouling from A-->C or from B-->D
            elif (sp2.Pcode[1] == "A" and sp1.Pcode[1] == "C") or (
                sp2.Pcode[1] == "B" and sp1.Pcode[1] == "D"
            ):
                process = "biofouling"
                if process in sp2.Pcompartment.processess:
                    sol = "k_" + process
                else:
                    sol=0

            # Defouling from C-->A or from D-->B
            elif (sp2.Pcode[1] == "C" and sp1.Pcode[1] == "A") or (
                sp2.Pcode[1] == "D" and sp1.Pcode[1] == "B"
            ):
                process = "defouling"
                if process in sp2.Pcompartment.processess:
                    sol = "k_" + process
                else:
                    sol=0
            else:
                sol=0
        else:
            sol=0

    # Different compartments--> Transport processess
    # settling, rising, mixing, resusp, advective transport, difussion, runoff, percolation?

    # if compartments are in the list of compartment connexions select process of connexion for compartment and assign rate constant, else process has rate of 0

    elif sp1.Pcompartment.Cname in sp2.Pcompartment.connexions:
        #transport between compartments only for same aggregation state and same particle size
        if sp1.Pcode[:2] == sp2.Pcode[:2]:
            process = sp2.Pcompartment.connexions[sp1.Pcompartment.Cname]
            if type(process) == list:
                sol2 = []
                for p in process:
                    sol2.append("k_" + p)
                sol = "+".join(sol2)
            else:
                sol = "k_" + process
        else:
            sol=0
    else:
        sol=0

    return sol


def interactionProcess_knames(sp1, interactions_df_knames, system_particle_object_list):
    sol = []
    for sp2 in system_particle_object_list:

        # Same particle in the same box and compartment (losses)
        if sp1.Pcode == sp2.Pcode:
            sol.append(interactions_df_knames[sp2.Pcode][sp1.Pcode])

        # Different particle or different river section or compartment
        else:

            # Same box (i.e. river section RS)--> In box processes

            if sp1.Pcode.split("_")[1] == sp2.Pcode.split("_")[1]:
                sol.append(inboxProcess_knames(sp1, sp2))

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

    # Asign loose rates
    elimination_rates = eliminationProcesses(system_particle_object_list, SpeciesList)

    interactions_df = pd.DataFrame(
        np.diag(elimination_rates), index=SpeciesList, columns=SpeciesList
    )

    # Asign interactions rates
    interactions_df_rows = []

    for sp1 in system_particle_object_list:
        interactions_df_rows.append(
            interactionProcess(sp1, interactions_df, system_particle_object_list)
        )

    # interact3(sp1) for sp1 in interactions_df.index.to_list()]
    array = np.column_stack(
        interactions_df_rows
    )  # vstack it was set as column stack and was wrong!!
    interactions_df_sol = pd.DataFrame(array, index=list_sp1, columns=list_sp1)

    return interactions_df_sol
