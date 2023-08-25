# This script contains the defined comartment objects for UTOPIA
# - List of COMPARTMENTS included:
#   air
#   oceanSurfaceWater
#   oceanMixedWater,
#   oceanColumnWater,
#   coastSurfaceWater,
#   coastColumnWater,
#   freshWaterSurface,
#   freshWaterBulk,
#   sediment,
#   urbanSoilSurface,
#   urbanSoil,
#   backgroundSoilSurface,
#   backgroundSoil,
#   agriculturalSoilSurface,
#   agriculturalSoil,

from objects.compartmetSubclasess import *


def UTOPIA_compartment_objects():

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

    oceanColumnWater = compartment_oceanWater(
        Cname="Ocean Column Water",
        SPM_mgL=1,
        waterFlow_m3_s=22100000,
        flowVelocity_m_s=0,
        T_K=278,
        G=1,
        Cvolume_m3=1e17,
        Cdepth_m=1000,
    )

    coastSurfaceWater = compartment_oceanWater(
        Cname="Coast Surface Water",
        SPM_mgL=7,
        T_K=278,
        G=10,
        flowVelocity_m_s=2.5,
        waterFlow_m3_s=0,
        Cvolume_m3=3.8e15,
        Cdepth_m=10,
    )

    coastColumnWater = compartment_oceanWater(
        Cname="Coast Column Water",
        SPM_mgL=7,
        T_K=278,
        G=1,
        waterFlow_m3_s=6840000000,
        flowVelocity_m_s=0,
        Cvolume_m3=1e17,
        Cdepth_m=10,
    )

    freshWaterSurface = compartment_water(
        Cname="Surface Freshwater",
        SPM_mgL=30,
        T_K=278,
        G=10,
        flowVelocity_m_s=1.3,
        waterFlow_m3_s=0,
        Cvolume_m3=100,
        Cdepth_m=10,
    )

    freshWaterBulk = compartment_water(
        Cname="Bulk Freshwater",
        SPM_mgL=30,
        T_K=278,
        G=10,
        flowVelocity_m_s=1.3,
        waterFlow_m3_s="nan",
        Cvolume_m3=100,
        Cdepth_m=30,
        Cwidth_m=50,
    )

    sediment = compartment_sediment(Cname="Sediment", Cvolume_m3=100, Cdepth_m=10)

    urbanSoilSurface = compartment_soil_surface(
        Cname="Urban Soil Surface",
        earthworm_density_in_m3=0,
        infiltration_capacity=0,
        Qrunoff_m3=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    urbanSoil = compartment_soil(
        Cname="Urban Soil",
        infiltration_capacity=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    backgroundSoilSurface = compartment_soil_surface(
        Cname="Background Soil Surface",
        earthworm_density_in_m3=0,
        infiltration_capacity=0,
        Qrunoff_m3=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    backgroundSoil = compartment_soil(
        Cname="Background Soil",
        infiltration_capacity=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    agriculturalSoilSurface = compartment_soil_surface(
        Cname="Agricultural Soil Surface",
        earthworm_density_in_m3=0,
        infiltration_capacity=0,
        Qrunoff_m3=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    agriculturalSoil = compartment_soil(
        Cname="Agricultural Soil",
        infiltration_capacity=0,
        Cvolume_m3=100,
        soilPore_waterVolume_m3=0,
        Cdepth_m=10,
    )

    air = compartment_air(
        Cname="Air",
        T_K=278,
        wind_speed_m_s=5,
        I_rainfall_mm=0,
        Cvolume_m3=100,
        Cdepth_m=10,
    )

    # Establish CONNEXIONS: only listed those compartments wich will recieve particles from the define compartment. i.e. the ocean surface water compartment transports particles to the ocean mix layer through settling and to air through sea spray resuspension
    oceanSurfaceWater.connexions = {
        "Ocean Mixed Water": "settling",
        "Air": "sea_spray_aerosol",
    }

    oceanMixedWater.connexions = {
        "Ocean Surface Water": "rising",
        "Ocean Column Water": "settling",
    }

    oceanColumnWater.connexions = {
        "Ocean Mixed Water": "rising",
        "Sediment": "settling",
    }

    coastSurfaceWater.connexions = {
        "Air": "sea_spray_aerosol",
        "Ocean Surface Water": "advective_transport",
        "Coast Column Water": "settling",
    }

    coastColumnWater.connexions = {
        "Coast Surface Water": "rising",
        "Ocean Mixed Water": "advective_transport",
        "Sediment": "settling",
    }

    freshWaterSurface.connexions = {
        "Coast Surface Water": "advective_transport",
        "Bulk Freshwater": "settling",
    }

    freshWaterBulk.connexions = {
        "Surface Freshwater": "rising",
        "Coast Column Water": "advective_transport",
        "Sediment": "settling",
    }

    sediment.connexions = {
        "Bulk Freshwater": "sediment_resuspension",
        "Coast Column Water": "sediment_resuspension",
        "Ocean Column Water": "sediment_resuspension",
    }

    urbanSoilSurface.connexions = {
        "Air": "soil_air_resuspension",
        "Urban Soil": ["percolation", "tillage"],
        "Surface Freshwater": "runoff_transport",
        "Coast Surface Water": "runoff_transport",
    }

    urbanSoil.connexions = {"Urban Soil Surface": "tillage"}

    backgroundSoilSurface.connexions = {
        "Air": "soil_air_resuspension",
        "Surface Freshwater": "runoff_transport",
        "Coast Surface Water": "runoff_transport",
        "Background Soil": ["percolation", "tillage"],
    }

    backgroundSoil.connexions = {"Background Soil Surface": ["percolation", "tillage"]}

    agriculturalSoilSurface.connexions = {
        "Air": "soil_air_resuspension",
        "Surface Freshwater": "runoff_transport",
        "Coast Surface Water": "runoff_transport",
        "Background Soil": "runoff_transport",
        "Agricultural Soil": ["percolation", "tillage"],
    }

    agriculturalSoil.connexions = {"Agricultural Soil Surface": "tillage"}

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

    return compartments
