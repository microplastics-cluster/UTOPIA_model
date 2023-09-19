import matplotlib.pyplot as plt
import numpy as np

def plot_bySize(results_dict,comp_name):
    #results_dict=Results_comp_organiced
    for agg in results_dict[comp_name]:
        y=results_dict[comp_name][agg]["number_of_particles"]
        #x=results_dict[comp_name][agg]["species"]
        x=["0.05","0.5","5","50","500"]
        plt.bar(x,y)
        plt.yscale("log")
        plt.title(agg +" in "+ comp_name)
        plt.ylabel("number of particles")
        plt.xlabel("Size bin (nm)")
    
   
            
      
        
            
        