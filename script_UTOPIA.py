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

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

# Generate objects
# Boxes

boxName = "Utopia"
UTOPIA = Box(boxName)
print(f"The model box {boxName} has been created")

modelBoxes = [UTOPIA]
# modelBoxes=instantiateBoxes_from_csv(boxFile)
boxNames_list = [b.Bname for b in modelBoxes]

# Compartmets
"""Call read imput file function for compartments"""

compartments = instantiate_compartments(inputs_path + "\inputs_compartments.csv")

# compartments = instantiate_compartments_from_csv(
#     inputs_path + "\inputs_compartments.csv"
# )

# Establish connexions between compartments defining their interaction mechanism: only listed those compartments wich will recieve particles from the define compartment. i.e. the ocean surface water compartment transports particles to the ocean mix layer through settling and to air through sea spray resuspension

set_interactions(compartments, connexions_path_file="compartmentsInteractions.csv")

# Assign modelling code to compartmanes
for c in range(len(compartments)):
    compartments[c].Ccode = c + 1

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

# Estimate volume of UTOPIA box by adding volumes of the compartments addedd
# UTOPIA.calc_Bvolume_m3() #currently volume of soil and air boxess are missing, to be added to csv file


# Add particles to compartments
for b in modelBoxes:
    for c in b.compartments:
        for p in particles:
            c.add_particles(copy.deepcopy(p))
    print(f"The particles have been added to the compartments of {b.Bname}")


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

fileName = "rateConstantsUTOPIA_Test.csv"

rate_constants_df = create_rateConstants_table(system_particle_object_list)

plot_rate_constants(rate_constants_df)

# "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed (k > 0.000556), while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.

RC_df_timeLim = timeLimit_RC(rate_constants_df, k=0.000556)

plot_rate_constants(RC_df_timeLim)

# rate_constants_df.loc[
#     rate_constants_df["k_heteroaggregation"] > 0.000556, "k_heteroaggregation"
# ] = 0.000556

# rate_constants_df.loc[
#     rate_constants_df["k_heteroaggregate_breackup"] > 0.000556,
#     "k_heteroaggregate_breackup",
# ] = 0.000556

# Save rate contants dataframe as csv file

df4 = RC_df_timeLim.fillna(0)
df4.to_csv(fileName, index=False)

# Plot rate constat values for comparison


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
