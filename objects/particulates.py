# Class Particulates generates particulate objects, especifically microplastic particle objects
# The class defines a particle object by its composition, shape and size

import math

from helpers.globalConstants import *


class Particulates:

    # constructor
    def __init__(
        self,
        Pname,
        Pform,
        Pcomposition,
        Pdensity_kg_m3,
        Pshape,
        PdimensionX_um,
        PdimensionY_um,
        PdimensionZ_um,
        t_half_d=5000,
        Pnumber=None,
    ):
        self.Pname = Pname
        self.Pform = Pform  # Pform has to be in the particles type list: ["freeMP",""heterMP","biofMP","heterBiofMP"]
        self.Pcomposition = Pcomposition
        self.Pdensity_kg_m3 = Pdensity_kg_m3
        self.Pshape = Pshape
        self.PdimensionX_um = PdimensionX_um  # shortest size
        self.PdimensionY_um = PdimensionY_um  # longest size
        self.PdimensionZ_um = PdimensionZ_um  # intermediate size
        self.PdimensionX_m = PdimensionX_um / 1000000  # shortest size
        self.PdimensionY_m = PdimensionY_um / 1000000  # longest size
        self.PdimensionZ_m = PdimensionZ_um / 1000000  # intermediate size
        self.Pnumber = Pnumber  # to be objetained from emissions and background concentration of the compartment
        self.radius_m = (
            self.PdimensionX_um / 1e6
        )  # In spherical particles from MP radius (x dimension)
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.Pemiss_t_y = 0  # set as 0 until the inputs_emissions.csv file is used
        self.t_half_d = t_half_d

    def __repr__(self):
        return (
            "{"
            + self.Pname
            + ", "
            + self.Pform
            + ", "
            + self.Pcomposition
            + ", "
            + self.Pshape
            + ", "
            + str(self.Pdensity_kg_m3)
            + ", "
            + str(self.radius_m)
            + "}"
        )

    # methods

    # volume calculation
    # different formulas for different particle shapes.
    # currently defined for spheres, fibres, cylinders, pellets and irregular fragments
    def calc_volume(self):

        if self.Pshape == "sphere":
            self.Pvolume_m3 = 4 / 3 * math.pi * (self.radius_m) ** 3
            # calculates volume (in m3) of spherical particles from MP radius (x dimension)
            self.CSF = 1
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        elif (
            self.Pshape == "fibre"
            or self.Pshape == "fiber"
            or self.Pshape == "cylinder"
        ):
            self.Pvolume_m3 = math.pi * (self.radius_m) ** 2 * (self.PdimensionY_m)
            # calculates volume (in m3) of fibres or cylinders from diameter and
            # length assuming cylindrical shape where X is the shorterst size (radius) ans Y the longest (heigth)
            self.CSF = (self.radius_m) / math.sqrt(self.PdimensionY_m * self.radius_m)
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        elif self.Pshape == "pellet" or self.Pshape == "fragment":
            self.Pvolume_m3 = (
                self.PdimensionX_m * self.PdimensionY_m * self.PdimensionZ_m
            )
            # approximate volume calculation for irregular fragments
            # approximated as a cuboid using longest, intermediate and shortest length
            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
            # might adjust later!!
            self.CSF = self.PdimensionX_m / math.sqrt(
                self.PdimensionY_m * self.PdimensionZ_m
            )
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        else:
            print("Error: unknown shape")
            # print error message for shapes other than spheres
            # (to be removed when other volume calculations are implemented)

    def calc_settling(self):
        # Assuming spherical shape (radius = X dimension)
        self.vSet_m_s = (
            2
            / 9
            * (self.density_kg_m3 - density_w_21C_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (self.radius_m) ** 2
        )

    def calc_numConc(self, concMass_mg_L, concNum_part_L):

        if concNum_part_L == 0:
            self.concNum_part_m3 = (
                concMass_mg_L / 1000 / self.Pdensity_kg_m3 / self.Pvolume_m3
            )
            # if mass concentration is given, it is converted to number concentration
        else:
            self.concNum_part_m3 = concNum_part_L * 1000
            # if number concentration is given, it is converted from part/L to part/m3

    def assign_compartment(self, comp):
        self.Pcompartment = comp

    # def calc_particleNum(self, Pemiss_t_y):

    #     if self.Pnumber == None:
    #         if self.Pemiss_t_y==
    #         self.Pnumber = Pemiss_t_y*1000/self.density_kg_m3/self.volume_m3
    #         #if mass concentration is given, it is converted to number concentration
    #     else:
    #         self.Pnumber = Pnumber
    #         #if number concentration is given, it is converted from part/L to part/m3
