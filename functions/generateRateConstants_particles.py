from functions import RC_generator


def generate_rateConstants(particle, spm, dict_comp):
    particle.RateConstants = dict.fromkeys(
        ["k_" + p for p in particle.Pcompartment.processess]
    )
    for process in particle.RateConstants:
        if (
            process[2:] == "heteroaggregation"
            or process[2:] == "heteroaggregate_breackup"
        ):
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, spm
            )

        elif process[2:] == "mixing":
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, dict_comp
            )

        elif process[2:] == "dry_depossition" or process[2:] == "wet_depossition":
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle, dict_comp
            )

        else:
            particle.RateConstants[process] = getattr(RC_generator, process[2:])(
                particle
            )
