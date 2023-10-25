# Class Compartment (parent class) generates compartment objects that belong by default to an assigned model box (Cbox)
# each compartment will contain four different particle objects corresponding to the 4 described aggregation states
import csv
import copy
from operator import attrgetter

MPforms_list = ["freeMP", "heterMP", "biofMP", "heterBiofMP"]


class Compartment:
    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        self.Cname = Cname
        self.Cdepth_m = Cdepth_m
        self.Clength_m = Clength_m
        self.Cwidth_m = Cwidth_m
        self.Cvolume_m3 = Cvolume_m3
        self.CsurfaceArea_m2 = CsurfaceArea_m2
        self.particles = {
            "freeMP": [],
            "heterMP": [],
            "biofMP": [],
            "heterBiofMP": [],
        }  # Each key corresponds to another dicttionary of size bins
        self.processess = [
            "degradation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
        ]
        self.connexions = []

    def assign_box(self, Box):
        self.CBox = Box

    def add_particles(self, particle):
        self.particles[particle.Pform].append(particle)
        particle.assign_compartment(self)

    def assign_particlesEmiss(self, emissionsFile):
        # Emissions are given as total in tons per year for each compartment but then distributed in percentages for the different size fractions
        # and MP types (default only emissions of pristine particles?)
        with open(emissionsFile, "r") as f:
            reader = csv.DictReader(f)
            emissions = list(reader)
        for i in emissions:
            if i["Cname"] == self.Cname:
                self.particles["freeMP"].sort(
                    key=attrgetter("radius_m")
                )  # Sorther from smallest to bigest radius. Migth need updating for other Pshapes
                self.particles["freeMP"][0].Pemiss_t_y = (
                    float(i.get("emissions_t_y")) * float(i.get("emiss_x01um_%")) / 100
                )
                self.particles["freeMP"][1].Pemiss_t_y = (
                    float(i.get("emissions_t_y")) * float(i.get("emiss_um_%")) / 100
                )
                self.particles["freeMP"][2].Pemiss_t_y = (
                    float(i.get("emissions_t_y")) * float(i.get("emiss_x10um_%")) / 100
                )
                self.particles["freeMP"][3].Pemiss_t_y = (
                    float(i.get("emissions_t_y")) * float(i.get("emiss_x100um_%")) / 100
                )
                self.particles["freeMP"][4].Pemiss_t_y = (
                    float(i.get("emissions_t_y")) * float(i.get("emiss_mm_%")) / 100
                )
            else:
                pass

    def assign_backgroundMPConc():
        pass  # Optionaly there can be background MPs concentrations of particle in the different compartments

    def calc_volume(self):
        if self.Cvolume_m3 is None:
            if any(
                attr is None for attr in [self.Cdepth_m, self.Clength_m, self.Cwidth_m]
            ):
                print(
                    "Missing parameters needded to calculate compartment volume --> Try calc_vol_fromBox or add missing values to compartment dimensions"
                )

            else:
                self.Cvolume_m3 = self.Cdepth_m * self.Clength_m * self.Cwidth_m
                # print(
                #     "Calculated "
                #     + self.Cname
                #     + " volume: "
                #     + str(self.Cvolume_m3)
                #     + " m3"
                # )
        else:
            pass
            # print("Assigned " + self.Cname + " volume: " + str(self.Cvolume_m3) + " m3")

    def calc_vol_fromBox(self):
        self.Cvolume_m3 = (
            self.CBox.Bvolume_m3 * self.CBox.CvolFractionBox[self.Cname.lower()]
        )

    def calc_particleConcentration_Nm3_initial(self):
        for p in self.particles:
            for s in self.particles[p]:
                self.particles[p][s].initial_conc_Nm3 = (
                    self.particles[p][s].Pnumber / self.Cvolume_m3
                )
