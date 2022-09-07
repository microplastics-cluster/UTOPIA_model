from objects.compartment import *

# Subclasses (inheritances) of the class compartment add extra attributes to the compatment that define the type of compartment


"""Would it be useful to have one class per specific compartment already parameterised to buil the model?
for one per compartment listed above for the UTOPIA model?"""

# UTOPIA_compartments = ["Ocean Surface Water","Ocean Mixed Water",
# "Ocean Column Water","Coast Surface Water","Coast Column Water","Surface Freshwater ",
# "Bulk FreshWater","Sediment","Urban Soil Surface","Urban Soil","Background Soil Surface",
# "Background Soil","Agricultural Soil Surface","Agricultural Soil","Air"]


# listOfProcessess=['k_discorporation','k_fragmentation','k_heteroaggregation',
# 'k_heteroaggregate_breackup','k_biofouling', 'k_defouling', 'k_advective_transport',
#  'k_settling', 'k_rising', "k_sea_spray_aerosol", "k_wind_trasport",
#  "k_wet_depossition", "k_dry_deposition", "k_sediment_resuspension",
#  "k_burial", "k_percolation", "k_runoff_transport", "k_tillage"]


class compartment_water(Compartment):
    # added new attributes relative to the water compartments such as SPM concentration, flow velocity etc.
    def __init__(
        self,
        Cname,
        SPM_mgL,
        flowVelocity_m_s,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):

        super().__init__(Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3)
        self.SPM_mgL = SPM_mgL
        self.flowVelocity_m_s = flowVelocity_m_s
        self.waterFlow_m3_s = waterFlow_m3_s
        self.T_K = T_K
        self.G = G  # Shear rate (G, in s−1)
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
        ]
        if flowVelocity_m_s == "nan":
            self.flowVelocity_m_s = self.waterFlow_m3_s / (
                self.Cdepth_m * self.Cwidth_m
            )


class compartment_oceanWater(compartment_water):
    # added new processess to the list of processess. new attributes that migth be needed to this processess should be added here
    def __init__(
        self,
        Cname,
        SPM_mgL,
        flowVelocity_m_s,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):
        super().__init__(
            Cname,
            SPM_mgL,
            flowVelocity_m_s,
            waterFlow_m3_s,
            T_K,
            G,
            Cdepth_m,
            Clength_m,
            Cwidth_m,
            Cvolume_m3,
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
            "sea_spray_aerosol",
        ]


class compartment_sediment(Compartment):
    def __init__(
        self, Cname, Cdepth_m=None, Clength_m=None, Cwidth_m=None, Cvolume_m3=None
    ):
        super().__init__(Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3)
        self.processess = [
            "discorporation",
            "fragmentation",
            "sediment_resuspension",
            "burial",
        ]


class compartment_soil(Compartment):
    def __init__(
        self,
        Cname,
        soilPore_waterVolume_m3,
        infiltration_capacity=0.25,  # from SimpleBox(4plastics)
        precipitation_rate=2.22 * 1**-8,  # from SimpleBox(4plastics)
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):
        super().__init__(Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3)
        self.processess = ["discorporation", "fragmentation", "percolation", "tillage"]
        self.infiltration_capacity = infiltration_capacity
        self.precipitation_rate = (precipitation_rate,)
        self.soilPore_waterVolume_m3 = soilPore_waterVolume_m3


class compartment_soil_surface(compartment_soil):
    def __init__(
        self,
        Cname,
        soilPore_waterVolume_m3,
        earthworm_density_in_m3,  # migth be needed for bioturbation
        Qrunoff_m3,  # migth be needed for runoff transport
        infiltration_capacity=0.25,
        precipitation_rate=2.22 * 1**-8,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):
        super().__init__(
            Cname,
            soilPore_waterVolume_m3,
            infiltration_capacity,
            precipitation_rate,
            Cdepth_m,
            Clength_m,
            Cwidth_m,
            Cvolume_m3,
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "runoff_transport",
            "tillage",
            "percolation",
            "soil_air_resuspension",
        ]
        self.earthworm_density_in_m3 = earthworm_density_in_m3
        self.Qrunoff_m3 = Qrunoff_m3


class compartment_deep_soil(compartment_soil):
    def __init__(
        self,
        Cname,
        soilPore_waterVolume_m3,
        infiltration_capacity=0.25,
        precipitation_rate=2.22 * 1**-8,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):
        super().__init__(
            Cname,
            soilPore_waterVolume_m3,
            infiltration_capacity,
            precipitation_rate,
            Cdepth_m,
            Clength_m,
            Cwidth_m,
            Cvolume_m3,
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "retention_in_soil",
            "sequestration_deep_soils",
        ]


# retention_in_soil (straining?) of the particles in soil following heteroaggregation with geocolloids?
# shall we also include heteroaggregation/heteroaggegrate break-up processess in the soil compartment? In SimpleBox for Nano they do account for aggregation and attachment

# Difference between retention in soil and sequestration deep soil: sequestrations deep soil is like burial in deep sediments (elemination process-->out of the system)


class compartment_air(Compartment):
    def __init__(
        self,
        Cname,
        T_K,
        wind_speed_m_s,
        I_rainfall_mm,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
    ):
        super().__init__(Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3)
        self.T_K = T_K
        self.wind_speed_m_s = wind_speed_m_s
        self.I_rainfall_mm = I_rainfall_mm
        self.processess = [
            "discorporation",
            "fragmentation",
            "wind_trasport",
            "dry_depossition",
            "wet_depossition",
        ]
        # shall we also include heteroaggregation/heteroaggegrate break-up processess in the air compartment?


class compartment_FullMulti_water(Compartment):
    def __init__(
        self, Cname, Cdepth_m=None, Clength_m=None, Cwidth_m=None, Cvolume_m3=None
    ):
        super().__init__(Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3)


# distinguish between dissolution and discorporation (agein?) and add as an extra process?
