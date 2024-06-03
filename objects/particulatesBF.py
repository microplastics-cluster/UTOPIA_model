from objects.particulates import Particulates  # class to generate MP and SPM object
from helpers import *


# define ParticulatesBF class (inheriting from Particulates class)
class ParticulatesBF(Particulates):
    "This is a class to create ParticulatesBIOFILM objects"

    # class attribute
    species = "particulate"

    # constructor
    def __init__(self, parentMP, spm):

        self.Pname = parentMP.Pname + "_BF"
        self.Pcomposition = parentMP.Pcomposition
        self.Pform = "biofMP"
        self.parentMP = parentMP
        self.BF_density_kg_m3 = spm.Pdensity_kg_m3
        self.BF_thickness_um = spm.PdimensionX_um
        self.radius_m = parentMP.radius_m + (
            self.BF_thickness_um / 1e6
        )  # In spherical particles from MP radius (x dimension)
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.t_half_d = 25000  # As per The Full Multi parameterization
        if parentMP.PdimensionY_um == 0:
            self.PdimensionY_um = 0
        else:
            self.PdimensionY_um = parentMP.PdimensionY_um + self.BF_thickness_um * 2

        if parentMP.PdimensionZ_um == 0:
            self.PdimensionZ_um = 0
        else:
            self.PdimensionZ_um = parentMP.PdimensionZ_um + self.BF_thickness_um * 2

        if parentMP.PdimensionX_um == 0:
            self.PdimensionX_um = 0
        else:
            self.PdimensionX_um = parentMP.PdimensionX_um + self.BF_thickness_um * 2

        self.Pshape = (
            parentMP.Pshape
        )  # to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)
        self.Pdensity_kg_m3 = (
            self.parentMP.radius_m**3 * self.parentMP.Pdensity_kg_m3
            + (
                (self.parentMP.radius_m + (self.BF_thickness_um / 1e6)) ** 3
                - self.parentMP.radius_m**3
            )
            * self.BF_density_kg_m3
        ) / ((self.parentMP.radius_m + (self.BF_thickness_um / 1e6)) ** 3)
        # equation from Kooi et al for density

        self.PdimensionX_m = self.PdimensionX_um / 1000000  # shortest size
        self.PdimensionY_m = self.PdimensionY_um / 1000000  # longest size
        self.PdimensionZ_m = self.PdimensionZ_um / 1000000  # intermediate size
