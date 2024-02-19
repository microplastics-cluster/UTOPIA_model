def overall_residence_time_calculation(tables_outputFlows, Results_extended):
    discorporation_flows = []
    for k in tables_outputFlows:
        discorporation_flows.append(sum(tables_outputFlows[k].k_discorporation))

    Tov_sec = sum(Results_extended["mass_g"]) / sum(discorporation_flows)

    Tov_days = Tov_sec / 86400
    Tov_years = Tov_days / 365

    print("Overall residence time (years): " + str(Tov_years))
