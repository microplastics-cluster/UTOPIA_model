"""Dynamic solver for UTOPIA"""
import numpy as np
import numpy.linalg
import pandas as pd
from scipy.integrate import solve_ivp
from helpers.helpers import *


def solve_ODES_Dynamic(
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
            PartMass_t0.at[sp_imput, "mass_g"] = imput_flows_g_s[sp_imput]

        # Input vector
        inputVector = PartMass_t0["mass_g"].to_list()

        # Get the rate constant matrix from the interactions DataFrame
        k = interactions_df.to_numpy()

        # Define the differential equation to solve to get the dynamic solution
        def dMdt(t, M, k, I):
            # To get a dynamic solution, we need to solve an ODE of the
            # form dM/dt = M.k + I, where M is the mass vector, k is the
            # constant matrix, and I is the emissions (input) vector
            # #TODO if we want inputs to be dynamic, we could make `I` a
            # function of time and pull out the relevant input value here
            dMdt = np.dot(M, k) + I
            return np.squeeze(dMdt)

        # Create time grid for 1 year, with units of seconds (because
        # inputs are in seconds - #TODO check this)
        # TODO also make the time span flexible
        t_eval = np.linspace(0, 86400 * 365, 366)
        t_span = (t_eval.min(), t_eval.max())
        inputs = inputVector

        # Use a daily timestep to solve, scaling up the inputVector
        # to g/day (from g/s)
        # t_eval = np.linspace(0, 365, 366)
        # t_span = (0, 365)
        # inputs = np.array(inputVector) * 86400

        res = solve_ivp(dMdt,
                        t_span=t_span,
                        y0=inputs,
                        args=(k, inputs),
                        t_eval=t_eval,
                        method='LSODA')
        print(res)
        np.save('temp_t', res.t)
        np.save('temp_y', res.y)

        # For the moment, just pull out the final timestep PECs
        # TODO add ability to save full timeseries (over and above
        # dumping the arrays to temporary files above)
        M_tfinal = res.y[:, -1]
        Results = pd.DataFrame({"species": SpeciesList, "mass_g": M_tfinal})

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
    ):

        print('ERROR: Particle number inputs not implemented for dynamic solutions yet')
    else:
        print("ERROR: No particles have been input to the system")

    return R, PartMass_t0
