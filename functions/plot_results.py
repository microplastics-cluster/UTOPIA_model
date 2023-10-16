import matplotlib.pyplot as plt
import numpy as np

particle_sizes_coding = {"mp5": "a", "mp4": "b", "mp3": "c", "mp2": "d", "mp1": "e"}

def plot_bySize_total_number_particles(results_dict,comp_name,dict_size_coding):
    
    new_size_dict=dict(zip([particle_sizes_coding[x] for x in dict_size_coding],[str(y) for y in dict_size_coding.values()]))
   
    fig, axs=plt.subplots(nrows=1,ncols=4,figsize=(15,7),sharey=True)

    for i, agg in enumerate(results_dict[comp_name]):
        y=results_dict[comp_name][agg]["number_of_particles"]
        X=results_dict[comp_name][agg]["species"]
        X_1=[s[0] for s in X]
        x=[new_size_dict[s] for s in X_1]
        axs[i].bar(x,y)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel("Total Number of Particles")
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)
    
   
def plot_bySize_total_mass(results_dict,comp_name,dict_size_coding):
    
    new_size_dict=dict(zip([particle_sizes_coding[x] for x in dict_size_coding],[str(y) for y in dict_size_coding.values()]))
   
    fig, axs=plt.subplots(nrows=1,ncols=4,figsize=(15,7),sharey=True)

    for i, agg in enumerate(results_dict[comp_name]):
        y=results_dict[comp_name][agg]["mass_g"]
        X=results_dict[comp_name][agg]["species"]
        X_1=[s[0] for s in X]
        x=[new_size_dict[s] for s in X_1]
        axs[i].bar(x,y)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel("Total mass (g)")
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)
      
        
            
        