import matplotlib.pyplot as plt


def generate_system_species_list(
    system_particle_object_list, MPforms_list, compartmentNames_list, boxNames_list
):
    particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}

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


import numpy as np
import pandas as pd


def timeLimit_RC(RC_df, k):
    processList = processList = [k for k in RC_df.columns if "k" in k]

    frame = np.asarray(RC_df[processList])
    frame[frame > k] = k

    RC_df = pd.DataFrame(frame, columns=processList)

    return RC_df
