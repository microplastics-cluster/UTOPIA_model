from functions import RC_generator


def generate_rateConstants(particle, spm, dict_comp, fsd, process_inputs_df):
    particle.RateConstants = dict.fromkeys(
        ["k_" + p for p in particle.Pcompartment.processess]
    )
    for process in particle.RateConstants:
        if (
            process[2:] == "heteroaggregation"
            or process[2:] == "heteroaggregate_breackup"
        ):
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, spm, process_inputs_df
            )

        elif process[2:] == "mixing":
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, dict_comp
            )

        elif process[2:] == "dry_deposition" or process[2:] == "wet_deposition":
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, dict_comp
            )
        elif process[2:] == "fragmentation":
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, fsd, process_inputs_df
            )
        elif (
            process[2:] == "discorporation"
            or process[2:] == "biofouling"
            or process[2:] == "defouling"
        ):
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, process_inputs_df
            )

        else:
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle
            )
