# UTOPIA documentation (to be translated to markdown)
## Model Structure
### OBJECTS
•	Model box (UTOPIA)
##### Attributes:
Bname   
Bdepth_m    
Blength_m   
Bwidth_m    
Bvolume_m3  
Bconexions  
Compartments =[]    

•	Model compartments

##### Shared compartment class attributes:

Cname   
Cdepth_m    
Clength_m   
Cwidth_m    
Cvolume_m3  
CsurfaceArea_m2 
particles = {"freeMP": [],"heterMP": [],"biofMP": [],"heterBiofMP":[]}Each key corresponds to another dicttionary of size bins    
connexions = [] 

UTOPIA_water_compartments = [
    "Ocean Surface Water",
    "Ocean Mixed Water",
    "Ocean Column Water",
    "Coast Surface Water",
    "Coast Column Water",
    "Surface Freshwater",
    "Bulk Freshwater",
]

##### Subclass specific attributes:
SPM_mgL 
waterFlow_m3_s  
T_K 
G   
flowVelocity_m_s    
processess = ["discorporation","fragmentation","heteroaggregation","heteroaggregate_breackup","biofouling","defouling","advective_transport","settling","rising","mixing"]  


UTOPIA_soil_compartments = [
    "Sediment",
    "Urban Soil Surface",
    "Urban Soil",
    "Background Soil Surface",
    "Background Soil",
    "Agricultural Soil Surface",
    "Agricultural Soil",
]

##### Subclass specific attributes:
infiltration_capacity=0.25,  # from SimpleBox(4plastics)    
precipitation_rate=2.22 * 1**-8,  # from SimpleBox(4plastics)   
soilPore_waterVolume_m3 
processess = ["discorporation","fragmentation","runoff_transport","tillage","percolation", "soil_air_resuspension"] 


UTOPIA_air_compartments = ["Air"]

##### Subclass specific attributes:
wind_speed_m_s  
I_rainfall_mm   
processess = ["discorporation","fragmentation","wind_trasport","dry_depossition","wet_depossition"] 

•	Particles

-	Particulates
-	ParticulatesBF
-	ParticulatesSPM

## MODEL MODIFICATIONS

•	COMPARTMENTS    

To add/change UTOPIA compartments go to compartments subclasses python file and add new compartment to the compartment list or water, soil or air compartment.
Add the compartment to the imput_compartments.csv file
If new attributes are added to the class this have to also be added in the class definition in the compartmentsSubclassess.py file as None attribute.

•	Rate constants

If we want to change rate constant values directly on the model these have to be changed at the particle level (not on the output table of Rate constants) so that the changes are reflected in the interactions matrix and model

Time limit can be included in the rate constants using the function timeLimit_particles_RC(system_particle_object_list, k), where k is the maximum value that k can take corresponding to it time limit 1/tlim (currently 30min on the processes that exceeds that speed (k > 0.000556)))


