# Simplified version of UTOPIA:

# - Only two particle SPECIES are defined:
#   Free plastic particles
#   Surface Modified particles (corresponding to particles that have been   aggregated with SPM and or biofouled)

# - List of COMPARTMENTS:
#   air
#   oceanSurfaceWater,
#   oceanColumnWater,
#   coastSurfaceWater,
#   coastColumnWater,
#   freshWaterSurface,
#   freshWaterBulk,
#   sediment,
#   backgroundSoilSurface,
#   backgroundSoil,
#   agriculturalSoilSurface,
#   agriculturalSoil

# -Run model for different size bins as series of consecutive steps (from bigger to smaller size bins)

import copy
import os
import pandas as pd
from functions import create_inputsTable_UTOPIA
from functions.create_rateConstants_tabel import *

from functions.create_inputsTable_UTOPIA import *

from functions.readImputs_from_csv import *
from helpers.helpers import *
from objects.box import *
from objects.UTOPIA_compartment_objects import *
from objects.compartmetSubclasess import *
from objects.particulates import *
from objects.particulatesBF import *
from objects.particulatesSPM import *

inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

## Generate objects

# MODEL BOX
"""Call read imput file function for model boxes (if more than one) in the same way as for MPs"""

boxName = "Utopia"
UTOPIA = Box(boxName)
print(f"The model box {boxName} has been created")

modelBoxes = [UTOPIA]
# modelBoxes=instantiateBoxes_from_csv(boxFile)
boxNames_list = [b.Bname for b in modelBoxes]

# COMPARTMENTS

compartments = UTOPIA_compartment_objects()


##Calculate compartments volume
for c in compartments:
    c.calc_volume()

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
