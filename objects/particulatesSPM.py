import math
from helpers import *
from objects.particulates import Particulates  # class to generate MP and SPM objects


# define ParticulatesSPM class (inheriting from Particulates class)
class ParticulatesSPM(Particulates):
    "This is a class to create ParticulatesSPM objects"

    # class attribute
    species = "particulate"

    # constructor
    def __init__(self, parentSPM, parentMP):

        self.Pname = parentMP.Pname + "_SPM"
        self.Pcomposition = parentMP.Pcomposition
        if parentMP.Pform == "biofMP":
            self.Pform = "heterBiofMP"
            self.t_half_d = 50000  # As per The Full multi parameterization
        else:
            self.Pform = "heterMP"
            self.t_half_d = 100000  # As per The Full multi parameterizatio
        self.parentMP = parentMP
        self.parentSPM = parentSPM
        self.Pdensity_kg_m3 = parentMP.Pdensity_kg_m3 * (
            parentMP.Pvolume_m3 / (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3)
        ) + parentSPM.Pdensity_kg_m3 * (
            parentSPM.Pvolume_m3 / (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3)
        )
        self.radius_m = (
            3 * (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3) / (4 * math.pi)
        ) ** (
            1 / 3
        )  # Note: this is an equivalent radius. MP-SPM most likely not truly spherical
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.Pshape = (
            parentMP.Pshape
        )  # to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)

    # methods

    # volume calculation - currently simple version.
    # more complexity to be added later:
    # different formulas for different particle shapes.
    # currently defined for spheres, fibres, cylinders, pellets and irregular fragments
    def calc_volume_heter(self, parentMP, parentSPM):
        if self.Pshape == "sphere":
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # calculates volume (in m3) of spherical particles from MP radius (x dimension)
            self.CSF = 1
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        elif (
            self.Pshape == "fibre"
            or self.Pshape == "fiber"
            or self.Pshape == "cylinder"
        ):
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # calculates volume (in m3) of fibres or cylinders from diameter and
            # length assuming cylindrical shape where X is the shorterst size (radius) ans Y the longest (heigth)
            self.CSF = (self.radius_m) / math.sqrt(self.PdimensionY_m * self.radius_m)
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        elif self.Pshape == "pellet" or self.Pshape == "fragment":
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # approximate volume calculation for irregular fragments
            # approximated as a cuboid using longest, intermediate and shortest length
            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
            # might adjust later!!
            self.CSF = self.PdimensionX_m / math.sqrt(
                self.PdimensionY_m * self.PdimensionZ_m
            )
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        else:
            print("Error: unknown shape")

        # print("Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3")
