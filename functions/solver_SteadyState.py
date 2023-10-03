import pandas as pd
from helpers.helpers import*

def solve_ODES_SS(system_particle_object_list,q_mass_g_s,q_num_s,sp_imput,interactions_df):
    SpeciesList=[p.Pcode for p in system_particle_object_list]
    
    if q_mass_g_s !=0:
        
        #set mass of particles for all particles in the system as zero
        m_t0=[]
        for p in system_particle_object_list:
            p.Pmass_g_t0=0
            m_t0.append(p.Pmass_g_t0)
        
        # dataframe of mass of particles at time 0
        PartMass_t0 = pd.DataFrame({"species": SpeciesList, "mass_g": m_t0})
        PartMass_t0 = PartMass_t0.set_index("species")
        
        # Set emissions
        PartMass_t0.at[sp_imput, "mass_g"] = -q_mass_g_s
        
        # Input vector
        inputVector = PartMass_t0["mass_g"].to_list()

        matrix = interactions_df.to_numpy()

        SteadyStateResults = np.linalg.solve(matrix, inputVector)

        Results = pd.DataFrame(
            {"species": SpeciesList, "mass_g": SteadyStateResults}
        )
        
        R=Results.set_index("species")
        for p in system_particle_object_list:
            p.Pmass_g_SS=R.loc[p.Pcode]["mass_g"]
            
        #Convert results in mass to particle number and add to the particle objects
        for p in system_particle_object_list:
            if "SPM" in p.Pname:
                if "BF" in p.Pname:
                    p.Pnum_SS=mass_to_num(mass_g=p.Pmass_g_SS,volume_m3=p.parentMP.parentMP.Pvolume_m3,density_kg_m3=p.parentMP.parentMP.Pdensity_kg_m3)
                else:
                    p.Pnum_SS=mass_to_num(mass_g=p.Pmass_g_SS,volume_m3=p.parentMP.Pvolume_m3,density_kg_m3=p.parentMP.Pdensity_kg_m3)    
            else:
                p.Pnum_SS=mass_to_num(mass_g=p.Pmass_g_SS,volume_m3=p.Pvolume_m3,density_kg_m3=p.Pdensity_kg_m3)        
                
        #Add to Results dataframe
        for p in system_particle_object_list:
            R.loc[p.Pcode,"number_of_particles"]=p.Pnum_SS  
    
    elif  q_num_s !=0:
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
        
        #Assign steady state (SS) results to paticles in particle number

        R=Results.set_index("species")
        for p in system_particle_object_list:
            p.Pnum_SS=R.loc[p.Pcode]["number_of_particles"]
        
        #Convert results in particle number to mass and add to the particle objects
        for p in system_particle_object_list:
            if "SPM" in p.Pname:
                if "BF" in p.Pname:
                    p.Pmass_SS_g=num_to_mass(number=p.Pnum_SS,volume_m3=p.parentMP.parentMP.Pvolume_m3,density_kg_m3=p.parentMP.parentMP.Pdensity_kg_m3)
                else:
                    p.Pmass_SS_g=num_to_mass(number=p.Pnum_SS,volume_m3=p.parentMP.Pvolume_m3,density_kg_m3=p.parentMP.Pdensity_kg_m3)  
            else:
                p.Pmass_SS_g=num_to_mass(number=p.Pnum_SS,volume_m3=p.Pvolume_m3,density_kg_m3=p.Pdensity_kg_m3)  
        
        #Add to Results dataframe
        for p in system_particle_object_list:
            R.loc[p.Pcode,"mass_g"]=p.Pmass_SS_g 
    
    return Results, R