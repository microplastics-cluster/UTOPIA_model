import pandas as pd
from helpers.helpers import *


def solve_ODES_SS(
    system_particle_object_list, q_num_s, imput_flows_g_s, interactions_df
):
    SpeciesList = [p.Pcode for p in system_particle_object_list]

    # Set initial mass of particles to 0
    # Set SS mass of particles to =???
    if sum(imput_flows_g_s.values()) != 0:
        # set mass of particles for all particles in the system as zero
        m_t0 = []
        for p in system_particle_object_list:
            p.Pmass_g_t0 = 0
            m_t0.append(p.Pmass_g_t0)

        # dataframe of mass of particles at time 0
        PartMass_t0 = pd.DataFrame({"species": SpeciesList, "mass_g": m_t0})
        PartMass_t0 = PartMass_t0.set_index("species")

        # Set emissions
        for sp_imput in imput_flows_g_s.keys():
            PartMass_t0.at[sp_imput, "mass_g"] = -imput_flows_g_s[sp_imput]

        # Input vector
        inputVector = PartMass_t0["mass_g"].to_list()

        matrix = interactions_df.to_numpy()

        SteadyStateResults = np.linalg.solve(matrix, inputVector)

        Results = pd.DataFrame({"species": SpeciesList, "mass_g": SteadyStateResults})

        R = Results.set_index("species")
        for p in system_particle_object_list:
            p.Pmass_g_SS = R.loc[p.Pcode]["mass_g"]

        # Convert results in mass to particle number and add to the particle objects
        for p in system_particle_object_list:
            if "SPM" in p.Pname:
                if "BF" in p.Pname:
                    p.Pnum_SS = mass_to_num(
                        mass_g=p.Pmass_g_SS,
                        volume_m3=p.parentMP.parentMP.Pvolume_m3,
                        density_kg_m3=p.parentMP.parentMP.Pdensity_kg_m3,
                    )
                else:
                    p.Pnum_SS = mass_to_num(
                        mass_g=p.Pmass_g_SS,
                        volume_m3=p.parentMP.Pvolume_m3,
                        density_kg_m3=p.parentMP.Pdensity_kg_m3,
                    )
            else:
                p.Pnum_SS = mass_to_num(
                    mass_g=p.Pmass_g_SS,
                    volume_m3=p.Pvolume_m3,
                    density_kg_m3=p.Pdensity_kg_m3,
                )

        # Add to Results dataframe
        for p in system_particle_object_list:
            R.loc[p.Pcode, "number_of_particles"] = p.Pnum_SS
        ### Estimate SS concentration and add to particles
        for p in system_particle_object_list:
            p.C_g_m3_SS = p.Pmass_g_SS / float(p.Pcompartment.Cvolume_m3)
            R.loc[p.Pcode, "concentration_g_m3"] = p.C_g_m3_SS
            p.C_num_m3_SS = p.Pnum_SS / float(p.Pcompartment.Cvolume_m3)
            R.loc[p.Pcode, "concentration_num_m3"] = p.C_num_m3_SS

    elif (
        q_num_s != 0
    ):  # By default the inputs are always given in mass this piece of code only needed if inputs given in particle numbers but this is has to be included in the imputs sections and updated to reflect the same structure as the mass inputs (list of inputs and not a single value)
        # Set number of particles for all particles in the system as zero
        N_t0 = []
        for p in system_particle_object_list:
            p.Pnumber = 0
            N_t0.append(p.Pnumber)

        # dataframe of number of particles at time 0
        PartNum_t0 = pd.DataFrame({"species": SpeciesList, "number_of_particles": N_t0})
        PartNum_t0 = PartNum_t0.set_index("species")

        # Set emissions
        PartNum_t0.at[sp_imput, "number_of_particles"] = -q_num_s

        # Input vector
        inputVector = PartNum_t0["number_of_particles"].to_list()
        matrix = interactions_df.to_numpy()

        SteadyStateResults = np.linalg.solve(matrix, inputVector)

        Results = pd.DataFrame(
            {"species": SpeciesList, "number_of_particles": SteadyStateResults}
        )

        # Assign steady state (SS) results to paticles in particle number

        R = Results.set_index("species")
        for p in system_particle_object_list:
            p.Pnum_SS = R.loc[p.Pcode]["number_of_particles"]

        # Convert results in particle number to mass and add to the particle objects
        for p in system_particle_object_list:
            if "SPM" in p.Pname:
                if "BF" in p.Pname:
                    p.Pmass_g_SS = num_to_mass(
                        number=p.Pnum_SS,
                        volume_m3=p.parentMP.parentMP.Pvolume_m3,
                        density_kg_m3=p.parentMP.parentMP.Pdensity_kg_m3,
                    )
                else:
                    p.Pmass_g_SS = num_to_mass(
                        number=p.Pnum_SS,
                        volume_m3=p.parentMP.Pvolume_m3,
                        density_kg_m3=p.parentMP.Pdensity_kg_m3,
                    )
            else:
                p.Pmass_g_SS = num_to_mass(
                    number=p.Pnum_SS,
                    volume_m3=p.Pvolume_m3,
                    density_kg_m3=p.Pdensity_kg_m3,
                )

        # Add to Results dataframe
        for p in system_particle_object_list:
            R.loc[p.Pcode, "mass_g"] = p.Pmass_g_SS

        ### Estimate SS concentration and add to particles
        for p in system_particle_object_list:
            p.C_g_m3_SS = p.Pmass_g_SS / float(p.Pcompartment.Cvolume_m3)
            R.loc[p.Pcode, "concentration_g_m3"] = p.C_g_m3_SS
            p.C_num_m3_SS = p.Pnum_SS / float(p.Pcompartment.Cvolume_m3)
            R.loc[p.Pcode, "concentration_num_m3"] = p.C_num_m3_SS

    else:
        print("ERROR: No particles have been input to the system")

    return R, PartMass_t0
