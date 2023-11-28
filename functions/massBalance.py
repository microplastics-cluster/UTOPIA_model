from functions.fillInteractions_df_fun_OOP import eliminationProcesses
from helpers.helpers import num_to_mass

import pandas as pd
import numpy as np


def massBalance(R, system_particle_object_list, q_mass_g_s):
    # Estimate looses: loss processess=[discorporation, burial]
    # Lossess also from fragmentation of the smallest size bin
    loss_processess = ["k_discorporation", "k_burial", "k_sequestration_deep_soils"]
    elimination_rates = []
    for p in system_particle_object_list:
        if p.Pcode[0] == "a":
            elimination_rates.append(
                sum(
                    [
                        p.RateConstants[e]
                        for e in loss_processess
                        if e in p.RateConstants
                    ]
                )
                + p.RateConstants["k_fragmentation"]
            )
        else:
            elimination_rates.append(
                sum(
                    [
                        p.RateConstants[e]
                        for e in loss_processess
                        if e in p.RateConstants
                    ]
                )
            )
    # mass at Steady state
    m_ss = R["mass_g"]

    # output flow
    out_flow_g_s = sum(elimination_rates * m_ss)

    print("Difference inflow-outflow = " + str(q_mass_g_s - out_flow_g_s))


def compartment_massBalance(
    comp,
    tables_outputFlows,
    PartMass_t0,
    comp_dict_inverse,
    dict_comp,
    tables_inputFlows,
):
    loss_processess = ["k_discorporation", "k_burial", "k_sequestration_deep_soils"]

    transfer_processes = [
        "k_advective_transport",
        "k_rising",
        "k_settling",
        "k_sea_spray_aerosol",
        "k_sediment_resuspension",
        "k_runoff_transport",
        "k_percolation",
        "k_tillage",
        "k_soil_air_resuspension",
        "k_wind_trasport",
        "k_dry_depossition",
        "k_wet_depossition",
        "k_mixing",
    ]
    comp_loss_processess = loss_processess + transfer_processes

    output_flows = tables_outputFlows[comp]
    output_flows_sum = output_flows.sum()

    out_flow_comp_g_s = sum(
        [
            val
            for proc, val in zip(output_flows_sum.index, output_flows_sum)
            if proc in comp_loss_processess
        ]
    )

    # input flow
    # Emissions
    for i, s in zip(PartMass_t0.index, PartMass_t0.values):
        if sum(s) != 0:
            if comp_dict_inverse[float(i[2:-7])] == comp:
                emiss_flow_g_s = -sum(s)
            else:
                emiss_flow_g_s = 0

    transport_input_flow = sum(tables_inputFlows[comp].sum())

    # Mass balance per compartment
    print(
        "Difference inflow-outflow in "
        + comp
        + " is = "
        + str(emiss_flow_g_s + transport_input_flow - out_flow_comp_g_s)
    )


def global_massBalance(q_mass_g_s, tables_outputFlows):
    # Estimate looses: loss processess=[discorporation, burial]
    # Lossess also from fragmentation of the smallest size bin
    loss_processess = ["k_discorporation", "k_burial", "k_sequestration_deep_soils"]
    output_flows = []
    for comp in tables_outputFlows:
        frag_flows = []
        for i in tables_outputFlows[comp].index:
            if i[0] == "a":
                frag_flows.append(tables_outputFlows[comp].loc[i, "k_fragmentation"])
            else:
                pass
        output_flows.append(
            sum(
                [
                    tables_outputFlows[comp][e]
                    for e in loss_processess
                    if e in tables_outputFlows[comp].columns
                ]
            ).sum()
            + sum(frag_flows)
        )

    print("Difference inflow-outflow = " + str(q_mass_g_s - sum(output_flows)))
