import copy
import os
import pandas as pd
from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *

from functions.create_inputsTable_UTOPIA import *

from functions.readImputs_from_csv import *
from helpers.helpers import *
from objects.box import *
from objects.compartmetSubclasess import *
from objects.particulates import *
from objects.particulatesBF import *
from objects.particulatesSPM import *

# Generate objects

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

# Boxes
"""Call read imput file function for model boxes (if more than one) in the same way as for MPs"""

boxName = "Utopia"
UTOPIA = Box(boxName)
print(f"The model box {boxName} has been created")

modelBoxes = [UTOPIA]
# modelBoxes=instantiateBoxes_from_csv(boxFile)
boxNames_list = [b.Bname for b in modelBoxes]

# Compartmets
"""Call read imput file function for compartments same way as for MPs"""

compartments = instantiate_compartments_from_csv(
    inputs_path + "\inputs_compartments.csv"
)


oceanSurfaceWater = compartment_oceanWater(
    Cname="Ocean Surface Water",
    SPM_mgL=1,
    waterFlow_m3_s=22100000,
    T_K=278,
    flowVelocity_m_s=0,
    G=1,
    Cvolume_m3=3.8e15,  # Value taken from SimpleBox4nano parameterization as volume of ocean water compartments surface layer
    Cdepth_m=10,
)
# Establish connexions: only listed those compartments wich will recieve particles from the define compartment. i.e. the ocean surface water compartment transports particles to the ocean mix layer through settling and to air through sea spray resuspension
oceanSurfaceWater.connexions = {
    "Ocean Mixed Water": "settling",
    "Air": "sea_spray_aerosol",
}

oceanMixedWater = compartment_oceanWater(
    Cname="Ocean Mixed Water",
    SPM_mgL=1,
    waterFlow_m3_s=22100000,
    flowVelocity_m_s=0,
    T_K=278,
    G=1,
    Cvolume_m3=1e16,
    Cdepth_m=100,
)
oceanMixedWater.connexions = {
    "Ocean Surface Water": "rising",
    "Ocean Column Water": "settling",
}

oceanColumnWater = compartment_oceanWater(
    Cname="Ocean Column Water",
    SPM_mgL=1,
    waterFlow_m3_s=22100000,
    flowVelocity_m_s=0,
    T_K=278,
    G=1,
    Cvolume_m3=1e18,
    Cdepth_m=1000,
)
oceanColumnWater.connexions = {"Ocean Mixed Water": "rising", "Sediment": "settling"}

coastSurfaceWater = compartment_oceanWater(
    Cname="Coast Surface Water",
    SPM_mgL=7,
    T_K=278,
    G=10,
    flowVelocity_m_s=2.5,
    waterFlow_m3_s=6840000000,
    Cvolume_m3=3.8e15,
    Cdepth_m=10,
)
coastSurfaceWater.connexions = {
    "Air": "sea_spray_aerosol",
    "Ocean Surface Water": "advective_transport",
    "Coast Column Water": "settling",
}


coastColumnWater = compartment_oceanWater(
    Cname="Coast Column Water",
    SPM_mgL=7,
    T_K=278,
    G=1,
    waterFlow_m3_s=6840000000,
    flowVelocity_m_s=2.5,
    Cvolume_m3=1e17,
    Cdepth_m=10,
)
coastColumnWater.connexions = {
    "Coast Surface Water": "rising",
    "Ocean Mixed Water": "advective_transport",
    "Sediment": "settling",
}

freshWaterSurface = compartment_water(
    Cname="Surface Freshwater",
    SPM_mgL=30,
    T_K=278,
    G=10,
    flowVelocity_m_s=1.3,
    waterFlow_m3_s=10000,  # placeholder
    Cvolume_m3=1.7e13,
    Cdepth_m=10,
)
freshWaterSurface.connexions = {
    "Coast Surface Water": "advective_transport",
    "Bulk Freshwater": "settling",
}

freshWaterBulk = compartment_water(
    Cname="Bulk Freshwater",
    SPM_mgL=30,
    T_K=278,
    G=10,
    flowVelocity_m_s=1.3,
    waterFlow_m3_s=10000,  # placeholder
    Cvolume_m3=1.6e14,
    Cdepth_m=30,
    Cwidth_m=50,
)
freshWaterBulk.connexions = {
    "Surface Freshwater": "rising",
    "Coast Column Water": "advective_transport",
    "Sediment": "settling",
}

sediment = compartment_sediment(Cname="Sediment", Cvolume_m3=100, Cdepth_m=10)
sediment.connexions = {
    "Bulk Freshwater": "sediment_resuspension",
    "Coast Column Water": "sediment_resuspension",
    "Ocean Column Water": "sediment_resuspension",
}

urbanSoilSurface = compartment_soil_surface(
    Cname="Urban Soil Surface",
    earthworm_density_in_m3=0,
    infiltration_capacity=0,
    Qrunoff_m3=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
urbanSoilSurface.connexions = {
    "Air": "soil_air_resuspension",
    "Urban Soil": ["percolation", "tillage"],
    "Surface Freshwater": "runoff_transport",
    "Coast Surface Water": "runoff_transport",
}

urbanSoil = compartment_soil(
    Cname="Urban Soil",
    infiltration_capacity=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
urbanSoil.connexions = {"Urban Soil Surface": "tillage"}

backgroundSoilSurface = compartment_soil_surface(
    Cname="Background Soil Surface",
    earthworm_density_in_m3=0,
    infiltration_capacity=0,
    Qrunoff_m3=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
backgroundSoilSurface.connexions = {
    "Air": "soil_air_resuspension",
    "Surface Freshwater": "runoff_transport",
    "Coast Surface Water": "runoff_transport",
    "Background Soil": ["percolation", "tillage"],
}

backgroundSoil = compartment_soil(
    Cname="Background Soil",
    infiltration_capacity=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
backgroundSoil.connexions = {"Background Soil Surface": ["percolation", "tillage"]}

agriculturalSoilSurface = compartment_soil_surface(
    Cname="Agricultural Soil Surface",
    earthworm_density_in_m3=0,
    infiltration_capacity=0,
    Qrunoff_m3=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
agriculturalSoilSurface.connexions = {
    "Air": "soil_air_resuspension",
    "Surface Freshwater": "runoff_transport",
    "Coast Surface Water": "runoff_transport",
    "Background Soil": "runoff_transport",
    "Agricultural Soil": ["percolation", "tillage"],
}

agriculturalSoil = compartment_soil(
    Cname="Agricultural Soil",
    infiltration_capacity=0,
    Cvolume_m3=100,
    soilPore_waterVolume_m3=0,
    Cdepth_m=10,
)
agriculturalSoil.connexions = {"Agricultural Soil Surface": "tillage"}

air = compartment_air(
    Cname="Air", T_K=278, wind_speed_m_s=5, I_rainfall_mm=0, Cvolume_m3=100, Cdepth_m=10
)
air.connexions = {
    "Agricultural Soil Surface": ["dry_depossition", "wet_depossition"],
    "Background Soil Surface": ["dry_depossition", "wet_depossition"],
    "Urban Soil Surface": ["dry_depossition", "wet_depossition"],
    "Surface Freshwater": ["dry_depossition", "wet_depossition"],
    "Coast Surface Water": ["dry_depossition", "wet_depossition"],
    "Ocean Surface Water": ["dry_depossition", "wet_depossition"],
}


compartments = [
    oceanSurfaceWater,
    oceanMixedWater,
    oceanColumnWater,
    coastSurfaceWater,
    coastColumnWater,
    freshWaterSurface,
    freshWaterBulk,
    sediment,
    urbanSoilSurface,
    urbanSoil,
    backgroundSoilSurface,
    backgroundSoil,
    agriculturalSoilSurface,
    agriculturalSoil,
    air,
]

# Assign modelling code to compartmanes
for c in range(len(compartments)):
    compartments[c].Ccode = c + 1


print(f"The compartments {[c.Cname for c in compartments]} have been generated")
##Calculate compartments volume
for c in compartments:
    c.calc_volume()

## Dictionary of compartments
dict_comp = {
    item.Cname: item for item in compartments
}  # Before the compartments association to RS...migth need to come after to also reflect the CBox connexion

compartmentNames_list = [item.Cname for item in compartments]

# PARTICLES
MPforms_list = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]
##Free microplastics (freeMP)
MP_freeParticles = instantiateParticles_from_csv(
    inputs_path + "\inputs_microplastics.csv"
)
print(
    f"The free MP particles {[p.Pname for p in MP_freeParticles]} have been generated"
)

###Calculate freeMP volume
for i in MP_freeParticles:
    i.calc_volume()
    # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")


##Biofouled microplastics (biofMP)
MP_biofouledParticles = []
for i in MP_freeParticles:
    MP_biofouledParticles.append(ParticulatesBF(parentMP=i))
print(
    f"The biofouled MP particles {[p.Pname for p in MP_biofouledParticles]} have been generated"
)

###Calculate biofMP volume
for i in MP_biofouledParticles:
    i.calc_volume()
    # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")

##Heteroaggregated microplastics (heterMP)
spm = Particulates(
    Pname="spm1",
    Pform="suspendedParticulates",
    Pcomposition="Mixed",
    Pdensity_kg_m3=2000,
    Pshape="sphere",
    PdimensionX_um=0.5 / 2,
    PdimensionY_um=0,
    PdimensionZ_um=0,
)
spm.calc_volume()
print(f"spm Volume: {spm.Pvolume_m3} m3")
print(f"Density of spm: {spm.Pdensity_kg_m3} kg_m3")

MP_heteroaggregatedParticles = []
for i in MP_freeParticles:
    MP_heteroaggregatedParticles.append(ParticulatesSPM(parentMP=i, parentSPM=spm))
print(
    f"The heteroaggregated MP particles {[p.Pname for p in MP_heteroaggregatedParticles]} have been generated"
)

###Calculate heterMP volume
for i in MP_heteroaggregatedParticles:
    i.calc_volume_heter(i.parentMP, spm)
    # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")

##Biofouled and Heteroaggregated microplastics (biofHeterMP)
MP_biofHeter = []
for i in MP_biofouledParticles:
    MP_biofHeter.append(ParticulatesSPM(parentMP=i, parentSPM=spm))
# for i in MP_biofHeter:
#     print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")
print(
    f"The biofouled and heteroaggregated MP particles {[p.Pname for p in MP_biofHeter]} have been generated"
)

###Calculate biofHeterMP volume
for i in MP_biofHeter:
    i.calc_volume_heter(i.parentMP, spm)

particles = (
    MP_freeParticles
    + MP_biofouledParticles
    + MP_heteroaggregatedParticles
    + MP_biofHeter
)

particles_properties = {
    "Particle": ([p.Pname for p in particles]),
    "Radius_m": ([p.radius_m for p in particles]),
    "Volume_m3": ([p.Pvolume_m3 for p in particles]),
    "Density_kg_m3": ([p.Pdensity_kg_m3 for p in particles]),
    "Corey Shape Factor": ([p.CSF for p in particles]),
}

particles_df = pd.DataFrame(data=particles_properties)
# print(particles_df)
particles_df.to_csv("Particles_properties_output.csv", index=False)


# Assign compartmets to UTOPIA

for comp in compartments:
    UTOPIA.add_compartment(copy.deepcopy(comp))  # Check if the use of copy is correct!!

print(
    f"The compartments {[comp.Cname for comp in UTOPIA.compartments]} have been assigned to {UTOPIA.Bname } model box"
)


# Check if compartments and boxes (river sections) are consistent in dimensions (i.e. sum of volumes of the different compartments can not exceed the volume of the box)

"""To be implemented for UTOPIA"""

# Add particles to compartments
for b in modelBoxes:
    for c in b.compartments:
        for p in particles:
            c.add_particles(copy.deepcopy(p))
    print(
        f"The particles have been added to the compartments of the river section {b.Bname}"
    )


# Based on the given model structure (created model boxes, compartments and particles)
# generate the process inputs table

process_inputs_df = create_inputsTable_UTOPIA(compartments, modelBoxes, inputs_path)

"""Revisit create inputs table function...assumptions to be discussed and parameters to be added"""


# Estimate rate constants per particle
from functions.generateRateConstants_particles import *

# List of particle objects in the system:
system_particle_object_list = []

for b in modelBoxes:
    for c in b.compartments:
        for freeMP in c.particles["freeMP"]:
            system_particle_object_list.append(freeMP)
        for heterMP in c.particles["heterMP"]:
            system_particle_object_list.append(heterMP)
        for biofMP in c.particles["biofMP"]:
            system_particle_object_list.append(biofMP)
        for heterBiofMP in c.particles["heterBiofMP"]:
            system_particle_object_list.append(heterBiofMP)

# Generate list of species names and add code name to object
SpeciesList = generate_system_species_list(
    system_particle_object_list, MPforms_list, compartmentNames_list, boxNames_list
)

for particle in system_particle_object_list:
    generate_rateConstants(particle, spm, dict_comp)

# create rate constants table:

fileName = "rateConstantsUTOPIA_Test230619.csv"

rate_constants_df = create_rateConstants_table(system_particle_object_list)

# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

rate_constants_df.loc[
    rate_constants_df["k_heteroaggregation"] > 0.000556, "k_heteroaggregation"
] = 0.000556

rate_constants_df.loc[
    rate_constants_df["k_heteroaggregate_breackup"] > 0.000556,
    "k_heteroaggregate_breackup",
] = 0.000556

# Save rate contants dataframe as csv file

df4 = rate_constants_df.fillna(0)
df4.to_csv(fileName, index=False)

# Plot rate constat values for comparison

processList = [
    "k_discorporation",
    "k_fragmentation",
    "k_heteroaggregation",
    "k_heteroaggregate_breackup",
    "k_settling",
    "k_rising",
    "k_advective_transport",
    "k_mixing",
    "k_biofouling",
    "k_sediment_resuspension",
    "k_burial",
    "k_defouling",
]

df_RC = df4[processList]

df_RC.plot(
    title="Rate constant values (s-1)",
    subplots=True,
    figsize=(10, 15),
    sharex=True,
    fontsize=12,
    stacked=True,
)

# Generate system of differentia equations (1-Matrix of interactions, 2-System of differential equations)

# Build Matrix of interactions
from functions.fillInteractions_df_fun_OOP import *

interactions_df = fillInteractions_fun_OOP(system_particle_object_list, SpeciesList)

# """SOLVE SYSTEM OF ODES"""


"""I have removed the BOX from the compartments option so that when you assign compartments to a box it automatically assignes the box to that compartment...but to be fixed with the copy function!!"""

# """SOLVE SYSTEM OF ODES"""

# # Initial number of particles in the system

# Set number of particles for all particles in the system as zero
N_t0 = []
for p in system_particle_object_list:
    p.Pnumber = 0
    N_t0.append(p.Pnumber)

# dataframe of number of particles at time 0
PartNum_t0 = pd.DataFrame({"species": SpeciesList, "number_of_particles": N_t0})
PartNum_t0 = PartNum_t0.set_index("species")

# Set values !=0

PartNum_t0.at["eA1_Utopia", "number_of_particles"] = 100


# Input vector

inputVector = PartNum_t0["number_of_particles"].to_list()

matrix = interactions_df.to_numpy()

SteadyStateResults = np.linalg.solve(matrix, inputVector)

Results = pd.DataFrame(
    {"species": SpeciesList, "number_of_particles": SteadyStateResults}
)

# Check the result is correct

np.allclose(np.dot(matrix, SteadyStateResults), inputVector)
##I change the sign of Steady state results because our ODES: 0=M*C+I --> C=-I*M^(-1)

# # Vector of volumes corresponding to the compartments of the river
# dilution_vol_m3 = volumesVector(Clist, compartments_prop)

# ConcFinal_num_m3 = pd.DataFrame(data=0, index=t_span, columns=Clist)
# for ind in range(len(NFinal_num)):
#     ConcFinal_num_m3.iloc[ind] = NFinal_num.iloc[ind]/dilution_vol_m3
