from functions.fillInteractions_df_fun_OOP import eliminationProcesses
from helpers.helpers import num_to_mass

import pandas as pd
import numpy as np

def calculate_losses_and_sources(system_particle_object_list,SpeciesList,inputVector,interactions_df):
    # Estimate looses
    elimination_rates = eliminationProcesses(system_particle_object_list, SpeciesList)

    losses_matrix=pd.DataFrame(
            np.diag(elimination_rates), index=SpeciesList, columns=SpeciesList
        )

    M=losses_matrix.to_numpy()

    Losses = np.linalg.solve(M, inputVector)

    Losses_result = pd.DataFrame(
        {"species": SpeciesList, "number_of_particles": Losses}
    )

    #Estimate sources

    source_M=interactions_df

    source_M.values[[np.arange(source_M.shape[0])]*2] = 0

    MS=source_M.to_numpy()

    Sources = np.linalg.solve(MS, inputVector)

    Sources_result = pd.DataFrame(
        {"species": SpeciesList, "number_of_particles": Sources}
    )
    
    
    
def massBalance(Results,system_particle_object_list, PartNum_t0):
    for p in system_particle_object_list:
        p.Pmass_SS_g=num_to_mass(number=p.Pnum_SS,volume_m3=p.Pvolume_m3,density_kg_m3=p.Pdensity_kg_m3)
        
    for p in system_particle_object_list:
        p.Pmass_t0_g=num_to_mass(number=p.Pnumber_t0,volume_m3=p.Pvolume_m3,density_kg_m3=p.Pdensity_kg_m3)
        
    mass_t0_g=[]
    mass_SS_g=[]
    for s in system_particle_object_list:
        mass_t0_g.append(s.Pmass_t0_g)
        mass_SS_g.append(s.Pmass_SS_g)
    
    
    

def test_solution_mass_balance(matrix,SteadyStateResults,inputVector,q_num_s):
    
    x=[]
    for i in range(len(matrix)):
        a=matrix[i]
        b=SteadyStateResults
        x.append(np.dot(a,b))
    print(sum(x))#Has to be equal to the inputs
        
    
    
        
    