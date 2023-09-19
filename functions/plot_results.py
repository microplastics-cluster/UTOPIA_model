import matplotlib.pyplot as plt
import numpy as np

def plot_bySize(results_dict,comp_name):
    #x=Results_comp_organiced[comp_name][agg]["species"]
    x=["0.05","0.5","5","50","500"]
    fig, axs=plt.subplots(nrows=1,ncols=4,figsize=(15,7))

    for i, agg in enumerate(results_dict[comp_name]):
        y=results_dict[comp_name][agg]["number_of_particles"]
        axs[i].bar(x,y)
        axs[i].title.set_text(agg)
        axs[i].set_yscale("log")
        axs[i].set_ylabel("number of particles")
        axs[i].set_xlabel("Size bin (nm)")
    fig.suptitle(comp_name)
    
   
            
      
        
            
        