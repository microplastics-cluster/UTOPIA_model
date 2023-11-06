###Modify rate constants by stablishing a time limit or chaging specific rate constant values using the change_RC_value function
# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

from functions.massBalance import*
from functions.create_rateConstants_tabel import *
from functions.solver_SteadyState import*
from functions.extract_results import*
from functions.plot_results import*
from functions.fillInteractions_df_fun_OOP import *


def timeLimit_particles_RC(particle_object_list, lim,SpeciesList,q_mass_g_s,q_num_s,sp_imput,model_lists,MPforms_list):
    for particle in particle_object_list:
        for k in particle.RateConstants:
            if (
                particle.RateConstants[k] is not None
                and particle.RateConstants[k] > lim
            ):
                particle.RateConstants[k] = lim
            else:
                pass
    
    RC_df_timeLim = create_rateConstants_table(particle_object_list)

    #Plot rate constants
    plot_rate_constants(RC_df_timeLim)



    # Build Matrix of interactions

    interactions_df = fillInteractions_fun_OOP(particle_object_list, SpeciesList)


    # SOLVE SYSTEM OF ODES

    Results, R = solve_ODES_SS(system_particle_object_list=particle_object_list,q_mass_g_s=q_mass_g_s,q_num_s=q_num_s,sp_imput=sp_imput,interactions_df=interactions_df)

    
    massBalance(R,particle_object_list, q_mass_g_s)

    Results_comp_dict=extract_by_comp(R.reset_index(),model_lists["compartmentNames_list"])
    Results_comp_organiced=extract_by_aggSt(Results_comp_dict,MPforms_list)

    #Plot results
    particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


    for comp in Results_comp_organiced:
        plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])            



def change_RC_value(particle_object_list,rc_name,rc_val,SpeciesList,q_mass_g_s,q_num_s,sp_imput,model_lists,MPforms_list):
    for p in particle_object_list:
        for k in p.RateConstants:
            if k== rc_name:
                p.RateConstants[k]=rc_val
            else:
                pass

    RC_df_timeLim = create_rateConstants_table(particle_object_list)

    #Plot rate constants
    plot_rate_constants(RC_df_timeLim)



    # Build Matrix of interactions

    interactions_df = fillInteractions_fun_OOP(particle_object_list, SpeciesList)


    # SOLVE SYSTEM OF ODES

    Results, R = solve_ODES_SS(system_particle_object_list=particle_object_list,q_mass_g_s=q_mass_g_s,q_num_s=q_num_s,sp_imput=sp_imput,interactions_df=interactions_df)

    
    massBalance(R,particle_object_list, q_mass_g_s)

    Results_comp_dict=extract_by_comp(R.reset_index(),model_lists["compartmentNames_list"])
    Results_comp_organiced=extract_by_aggSt(Results_comp_dict,MPforms_list)

    #Plot results
    particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}


    for comp in Results_comp_organiced:
        plot_bySize_total_number_particles(Results_comp_organiced,comp,model_lists["dict_size_coding"])