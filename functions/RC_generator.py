# Extension of RS_generator module from the FUll Multi containing functions to calculate all rate constants


###So far programmed so that if a parameter value for estimating a rate constant is missing (i.e. NOne or nan), the rate constant is equal to zero....this has to be fixed!! Example,when no depth of a compartment is given the rate constant of burial or resuspension is equal to 0

import math
import pandas as pd
import os
import numpy as np


# import file storing required constants
from helpers.globalConstants import *

# # Read input data file

# process_inputs_df = pd.read_csv(
#     filepath_or_buffer=os.path.join(
#         os.path.dirname(__file__), "../inputs/processInputs_table.csv"
#     )
# )


def discorporation(particle, process_inputs_df):
    # degradation estimations
    # discorporation state will be given as output after runing
    # the model with no discorporation. degradation state will be given in time units as residence time in the compartment

    # Change process name from degradation
    # to discorporation from corporeal

    ##discorporation should be particle size dependent so that the higuer the surface area to volume ratio of the particel the fastest the process should happen (therfore for smaller particles fagmentation should happen at faster rates)

    # If we use the data from Pfohl et al. 2022 we would use a degradation rate of 6.3 x 10-6 but this is for particles of TPU-ether arom in the size range between 50-200um. We asume this value as discorporation rate for the 50 um MP plastics in free form

    # k_desint_50um_h= 6.3*10**-6 #(in h-1)
    # k_deg_free=k_desint_50um_h*(50/particle.diameter_um)

    """relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    # degradation half-life of MPs used as input is in days
    cond = (
        (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
    )
    t_half_d = float(process_inputs_df.loc[cond, "thalf_deg_d"])

    # degradation rate constant
    k_deg = math.log(2) / (t_half_d * 24 * 60 * 60)

    # * (
    #     50**2 / (particle.diameter_um) ** 2
    # )  # the discorporation rate is normalised to the surface area to volume ratio of the 50um particles since using t_half degradation rates derived from Pfohl et al. 2022 paper. this has already been included when building the table of input parameters

    return k_deg


def fragmentation(particle, fsd, process_inputs_df):

    # modelled as a size-dependent process based on an estimated rate constant (𝑘frag_gen= 1/tfrag_gen_d)
    # for fragmentation of pristine particles in the largest (x=5000μm => mp5 => e) size class.

    # estimate fragmentation relation between size bins using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html)

    # Fragmentation of heteroaggregated particles is assumed negligible in the default model formulation

    cond = (
        (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == "mp5")
    )
    t_frag_d = process_inputs_df.loc[cond, "tfrag_gen_d"].item()

    if t_frag_d == "NAN":
        # It is assumed that fragmentation doesnt occur in heteroaggregated particles (heterMP or heterBiofMP)
        frag_rate = 0
        fragments_formed = 0
    else:
        if (
            particle.Pname[0:3] == "mp1"
        ):  # print("Smallest sizeBin mp5(0.05um), fragments formed will be considered losses")
            frag_rate = 0  # We consider only discorporation for mp5(0.05um)
            fragments_formed = 0
        else:
            volume_fragment = (
                4 / 3 * math.pi * (float(particle.radius_m) / 10) ** 3
            )  #!!!only works for bins 10 times smaller!!!
            fragments_formed = float(particle.Pvolume_m3) / volume_fragment
            frag_rate = (
                (1 / (float(t_frag_d) * 24 * 60 * 60))
                * float(particle.diameter_um)
                / 1000
            )
    # each particle fractions into fragments of smaller sizes and the distribution is expresses via the fragment size distribution matrix fsd. # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class
    size_dict = {chr(i): i - ord("a") for i in range(ord("a"), ord("e") + 1)}

    k_frag = frag_rate * fsd[size_dict[particle.Pcode[0]]]

    return (
        k_frag.tolist()
    )  # I have removed the fragments formed from the output to have an homogeneus solution in the table of rate constants (consider dumping this values in another way later when/if needed (for Mass Balance?))


def settling(particle):
    # settling calculations
    """settling can be calculated using different equations (e.g. Stokes,
    modified versions of it or others) or be taken from experimental studies
    !! currently only classical Stokes is implemented (which is probably not
    very realistic and will be updated soon !!"""

    # Depending on the compartment we should use a specific water density

    if "Freshwater" in particle.Pcompartment.Cname:
        w_den_kg_m3 = density_w_21C_kg_m3
    else:
        w_den_kg_m3 = density_seaWater_kg_m3

    settlingMethod = "Stokes"

    # Settling occurs in all aquatic compartments which should be specified in the comprtment class
    # if particle.Pcompartment.Cname in ["Sediment", "Agricultural Soil","Urban Soil"...]
    #     k_set = 0

    if settlingMethod == "Stokes":
        vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - w_den_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (float(particle.radius_m)) ** 2
        )
    else:
        print("Error: cannot calculate settling other than Stokes yet")
        # print error message settling methods other than Stokes
        # (to be removed when other settling calculations are implemented)

    # for the water and surface water compartments:
    # settling and rising rate constants for free MP
    if vSet_m_s > 0:
        k_set = vSet_m_s / float(particle.Pcompartment.Cdepth_m)

    elif vSet_m_s < 0:
        k_set = 0

    else:
        k_set = 0

    return k_set


def rising(particle):
    # rising calculations
    """rising is calculated in the same way as settling for particles with negative
    settling velocities. It can be calculated using different equations (e.g. Stokes,
    modified versions of it or others) or be taken from experimental studies
    !! currently only classical Stokes is implemented (which is probably not
    very realistic and will be updated soon !!"""

    settlingMethod = "Stokes"

    # Rising only occus in the lower water compartments wich for UTOPIA are: ["Ocean Mixed Water",
    # "Ocean Column Water","Coast Column Water","Bulk FreshWater"]

    if particle.Pcompartment.Cname in [
        "Ocean_Mixed_Water",
        "Ocean_Column_Water",
        "Coast_Column_Water",
        "Bulk_Freshwater",
    ]:

        if "Freshwater" in particle.Pcompartment.Cname:
            w_den_kg_m3 = density_w_21C_kg_m3
        else:
            w_den_kg_m3 = density_seaWater_kg_m3

        if settlingMethod == "Stokes":
            vSet_m_s = (
                2
                / 9
                * (float(particle.Pdensity_kg_m3) - w_den_kg_m3)
                / mu_w_21C_kg_ms
                * g_m_s2
                * (float(particle.radius_m)) ** 2
            )
        else:
            print("Error: cannot calculate settling other than Stokes yet")
        # print error message settling methods other than Stokes
        # (to be removed when other settling calculations are implemented)
    else:
        vSet_m_s = 0
    # for the water and surface water compartments:
    # settling and rising rate constants for free MP
    if vSet_m_s > 0:
        k_rise = 0

    elif vSet_m_s < 0:
        k_rise = -vSet_m_s / float(particle.Pcompartment.Cdepth_m)

    else:
        k_rise = 0

    return k_rise


def heteroaggregation(particle, spm, process_inputs_df):
    if (particle.Pform == "freeMP") or (particle.Pform == "biofMP"):
        # heteroaggregation rate constants
        """heteroaggregation requires to particles to collide and interact
        favorably for the collision to result in attachment
        the heteroaggregation rate constants is therefore composed of two
        parts, 1) a collision rate constant and 2) and attachement
        efficiency (alpha) (representing the probability of attachement).
        For heteroaggregation a common simplifaction is the assumption that
        SPM concentration is not signficantly affected by the heteroaggre-
        gation process. Therefore, a pseudo first-order heteroaggregation
        rate constant is obtained by multiplying collision rate with alpha
        and with the SPM number concentration"""

        # first the different collision mechanisms are calculated
        k_peri = (
            (2 * k_B_J_K * float(particle.Pcompartment.T_K))
            / (3 * mu_w_21C_kg_ms)
            * (float(particle.radius_m) + spm.radius_m) ** 2
            / (float(particle.radius_m) * spm.radius_m)
        )
        # perikinetic contributions to collision rate constant (Brownian motion)

        k_ortho = (
            4
            / 3
            * float(particle.Pcompartment.G)
            * (float(particle.radius_m) + spm.radius_m) ** 3
        )
        # orthokinetic contributions to collision rate constant (caused by fluid motion)

        if "Freshwater" in particle.Pcompartment.Cname:
            w_den_kg_m3 = density_w_21C_kg_m3
        else:
            w_den_kg_m3 = density_seaWater_kg_m3

        MP_vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - w_den_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (float(particle.radius_m)) ** 2
        )

        SPM_vSet_m_s = (
            2
            / 9
            * (spm.Pdensity_kg_m3 - w_den_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (spm.radius_m) ** 2
        )
        # settling velocity. currently according to classical Stokes law. Need to include other modes and put calculation on its own, so that it can also be accessed for other processes

        k_diffSettling = (
            math.pi
            * (float(particle.radius_m) + spm.radius_m) ** 2
            * abs(MP_vSet_m_s - SPM_vSet_m_s)
        )

        # differential settling contributions to collision rate constant

        k_coll = k_peri + k_ortho + k_diffSettling
        # the collision rate constant

        cond_alpha = (
            (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
            & (process_inputs_df["MPform"] == particle.Pform)
            & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
        )
        alpha = process_inputs_df.loc[cond_alpha, "alpha_heter"].item()
        if alpha == "NAN":
            k_hetAgg = 0
        else:
            spm.calc_numConc(
                concMass_mg_L=float(particle.Pcompartment.SPM_mgL), concNum_part_L=0
            )
            SPM_concNum_part_m3 = spm.concNum_part_m3
            k_hetAgg = float(alpha) * k_coll * SPM_concNum_part_m3
            # the pseudo first-order heteroaggregation rate constant
    else:
        k_hetAgg = 0

    return k_hetAgg


def heteroaggregate_breackup(particle, spm, process_inputs_df):
    """Assumption: the breack-up of heteroaggregates is 10E8 times slower than the formation of heteroaggregates"""

    if (particle.Pform == "heterMP") or (particle.Pform == "heterBiofMP"):
        # Kbreackup is calculated based on Kheter of the free and biofouled MPs

        # data is limited on aggregate breakup, but this process is likely
        # more relvant for larger aggregates
        #!! 1/10 of k_hetAgg is just a placeholder,  needs to be refined
        # possibly using a size dependent function !!

        # first the different collision mechanisms are calculated

        k_peri = (
            (2 * k_B_J_K * float(particle.Pcompartment.T_K))
            / (3 * mu_w_21C_kg_ms)
            * (float(particle.radius_m) + spm.radius_m) ** 2
            / (float(particle.radius_m) * spm.radius_m)
        )
        # perikinetic contributions to collision rate constant (Brownian motion)

        k_ortho = (
            4
            / 3
            * float(particle.Pcompartment.G)
            * (float(particle.radius_m) + spm.radius_m) ** 3
        )
        # orthokinetic contributions to collision rate constant (caused by fluid motion)
        if "Freshwater" in particle.Pcompartment.Cname:
            w_den_kg_m3 = density_w_21C_kg_m3
        else:
            w_den_kg_m3 = density_seaWater_kg_m3

        MP_vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - w_den_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (float(particle.radius_m)) ** 2
        )
        SPM_vSet_m_s = (
            2
            / 9
            * (spm.Pdensity_kg_m3 - w_den_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (spm.radius_m) ** 2
        )
        # settling velocity. currently according to classical Stokes law. Need to include other modes and put calculation on its own, so that it can also be accessed for other processes

        k_diffSettling = (
            math.pi
            * (float(particle.radius_m) + spm.radius_m) ** 2
            * abs(MP_vSet_m_s - SPM_vSet_m_s)
        )
        # differential settling contributions to collision rate constant

        k_coll = k_peri + k_ortho + k_diffSettling
        # the collision rate constant
        if particle.Pform == "heterMP":
            cond_alpha = (
                (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
                & (process_inputs_df["MPform"] == "freeMP")
                & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
            )
        elif particle.Pform == "heterBiofMP":
            cond_alpha = (
                (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
                & (process_inputs_df["MPform"] == "biofMP")
                & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
            )

        alpha = process_inputs_df.loc[cond_alpha, "alpha_heter"].item()

        if alpha == "NAN":
            k_aggBreakup = 0
        else:
            spm.calc_numConc(
                concMass_mg_L=float(particle.Pcompartment.SPM_mgL), concNum_part_L=0
            )
            SPM_concNum_part_m3 = spm.concNum_part_m3
            k_hetAgg = float(alpha) * k_coll * SPM_concNum_part_m3
            # the pseudo first-order heteroaggregation rate constant

            k_aggBreakup = (1 / 1000000000) * k_hetAgg
    else:
        k_aggBreakup = 0

    return k_aggBreakup


def advective_transport(particle):
    k_adv = float(particle.Pcompartment.waterFlow_m3_s) / float(
        particle.Pcompartment.Cvolume_m3
    )

    return k_adv


def mixing(particle, dict_comp):

    # Now adapted to UTOPIA's compartments and changed rates
    # k_mix has to be multiplied by the compartment volume ratio calculated with the interacting compartment volume

    # k_mix_up = (
    #     10**-2
    # )  # (1): <Handbook of Chemical Mass Transport in the Environment> Edited by Louis J. Thibodeaux, Donald Mackay (DOI: 10.1201/b10262)

    # k_mix_down = (
    #     10**-3
    # )  # (2): <Handbook on Mixing in Rivers> Edited by J.C. Rutherford (Water and Soil Miscellaneous Publication No. 26. 1981. 60pp.ISSN 0110-4705)

    # Assuming that vertical mixing for the surface of the ocean and the coast is in the 1 hour time scale. The water in the surface will take 60 min to travel 50 m, half way trhough the mix layer (100m deep). 50m/60min = 0.83 m/min= 0.0138 m/s

    # FROM OECD tool: ocean water mixing to the depths of hell (Implies residence time in the mixed layer of 100 years, Wania and Mackay GloboPOP Value, Sci Tot Env. 1995) : 1.14 * 10 ^ -4 m/h=3.167E-8 m/s

    # Assuming that vertical mixing of freshwater compartments is in the minutes timescale. The water in the surface will take 2 min to travel 5 m, half way trhough the mix layer (10m deep). 2m/5min = 0.4 m/min= 0.0067 m/s

    flowRate_mixUP_ocean_m3_s = 0.0138 * float(
        dict_comp["Ocean_Mixed_Water"].CsurfaceArea_m2
    )

    flowRate_mixDown_ocean_m3_s = 3.167e-8 * float(
        dict_comp["Ocean_Mixed_Water"].CsurfaceArea_m2
    )

    flowRate_mix_coast_m3_s = 0.0138 * float(
        dict_comp["Coast_Column_Water"].CsurfaceArea_m2
    )

    flowRateMix_freshWater_m3_s = 0.0067 * float(
        dict_comp["Bulk_Freshwater"].CsurfaceArea_m2
    )

    if particle.Pcompartment.Cname == "Ocean_Mixed_Water":
        k_mix_up = flowRate_mixUP_ocean_m3_s / float(particle.Pcompartment.Cvolume_m3)
        k_mix_down = flowRate_mixDown_ocean_m3_s / float(
            particle.Pcompartment.Cvolume_m3
        )

        k_mix = [k_mix_up, k_mix_down]
        # {"mix_up": k_mix_up, "mix_down": k_mix_down}

    elif particle.Pcompartment.Cname == "Ocean_Column_Water":
        k_mix = flowRate_mixDown_ocean_m3_s / float(particle.Pcompartment.Cvolume_m3)

    elif particle.Pcompartment.Cname == "Ocean_Surface_Water":
        k_mix = flowRate_mixUP_ocean_m3_s / float(particle.Pcompartment.Cvolume_m3)

    elif particle.Pcompartment.Cname in ["Coast_Column_Water", "Coast_Surface_Water"]:
        k_mix = flowRate_mix_coast_m3_s / float(particle.Pcompartment.Cvolume_m3)
    elif particle.Pcompartment.Cname in ["Surface_Freshwater", "Bulk_Freshwater"]:
        k_mix = flowRateMix_freshWater_m3_s / float(particle.Pcompartment.Cvolume_m3)

    else:
        print("No mixing implemented for this compartment")
        k_mix = 0

    return k_mix


def biofouling(particle, process_inputs_df):
    cond_biof = (
        (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
    )
    t_biof_growth_d = process_inputs_df.loc[cond_biof, "tbiof_growth_d"].item()

    if t_biof_growth_d == "NAN":
        k_biof = 0
    else:
        k_biof = 1 / float(t_biof_growth_d) / 24 / 60 / 60

    # assume it takes x days for biofilm coverage to grow

    return k_biof


def defouling(particle, process_inputs_df):
    # Defouling = degradation of Biofilm.

    cond_defoul = (
        (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
    )
    tbiof_degrade_d = process_inputs_df.loc[cond_defoul, "tbiof_degrade_d"].item()
    if tbiof_degrade_d == "NAN":
        k_defoul = 0
    # assume it takes x days for biofilm coverage to be degraded
    else:
        k_defoul = 1 / float(tbiof_degrade_d) / 24 / 60 / 60

    return k_defoul


def sediment_resuspension(particle):
    # When no depth parameter available assign transfer sediment to water rate taken from SimpleBox for Plastics model
    # Currently placeholder values. To be revisited
    resusp_dict = {
        "Sediment_Freshwater": 1e-9,
        "Sediment_Coast": 1e-10,
        "Sediment_Ocean": 1e-11,
    }

    k_resusp = resusp_dict[particle.Pcompartment.Cname]

    return k_resusp


def burial(particle):
    # Currenlty place holder values. To be revisited

    # When no depth parameter available assign burail rate taken from SimpleBox for Plastics model
    burial_dict = {
        "Sediment_Freshwater": 2.7e-9,
        "Sediment_Coast": 1e-9,
        "Sediment_Ocean": 5e-10,
    }

    k_burial = burial_dict[particle.Pcompartment.Cname]

    return k_burial


def soil_air_resuspension(particle):
    # REF: BETR global approach for soil-air resuspension (10^-10 m/h)

    sar_rate = 10e-10 / 60 / 60

    sar_rate_dict = {
        "a": sar_rate,
        "b": sar_rate,
        "c": sar_rate,
        "d": sar_rate / 100,
        "e": sar_rate / 1e4,
    }

    k_sa_reusp = sar_rate / float(particle.Pcompartment.Cdepth_m)

    return k_sa_reusp


def soil_convection(particle):
    # Mixing of soil particles via bioturbation and freeze/thaw cycles

    MTCsconv = 4.54e-7
    # From the OECD Tool: MTCsconv = 4.54 * 10 ^-7 (m/h)'soil side solid phase convection MTC

    k_soil_convection = (MTCsconv / (60 * 60)) / float(particle.Pcompartment.Cdepth_m)

    # if particle.Pcompartment.Cname in [
    #     "Urban_Soil_Surface",
    #     "Background_Soil_Surface",
    #     "Agricultural_Soil_Surface",
    # ]:
    #     k_soil_convection = (
    #         (C_massTransfer_m_h /(60 * 60 ))/ float(particle.Pcompartment.Cdepth_m)
    #     )
    # elif particle.Pcompartment.Cname in [
    #     "Beaches_Deep_Soil",
    #     "Impacted_Soil",
    #     "Background_Soil",
    # ]:
    #     k_soil_conv = (C_massTransfer_m_h /( 60 * 60)) / float(particle.Pcompartment.Cdepth_m)
    #     k_soil_conv_down = (
    #         (1 / 20) * (C_massTransfer_m_h /( 60 * 60 )/ float(particle.Pcompartment.Cdepth_m))
    #     )
    #     k_soil_convection = [k_soil_conv, k_soil_conv_down]

    # else:
    #     k_soil_convection = 0

    return k_soil_convection


def percolation(particle):
    # downwards movement of particles in soil via infiltrated water
    # # to be defined/formulated

    # k_percol = particle.Pcompartment.infiltration_capacity*particle.Pcompartment.precipitation_rate*(float(particle.Pcompartment.Cvolume_m3)/float(particle.Pcompartment.Cdepth_m))/float(particle.Pcompartment.soilPore_waterVolume_m3)

    k_percol = 0

    return k_percol


def runoff_transport(particle):
    # transport from top soil layers to surface waters ["Coast_Surface_Water","Surface_Freshwater"] via runoff water
    # to be formulated

    # REF: BETR global approach for MTCsoilrunoff = 2.3 * 10 ^ -8  (m/h) 'soil solids runoff rate  (Scheringer, P230)

    runooff_dict = {
        "Beaches_Soil_Surface": 2.3e-8,
        "Background_Soil_Surface": 2.3e-8,
        "Impacted_Soil_Surface": 2.3e-8,
    }
    runoff_rate = (
        runooff_dict[particle.Pcompartment.Cname]
        / float(particle.Pcompartment.Cdepth_m)
    ) / (60 * 60)

    # The total amount of runoff will be distributed into the recieving compartments according to the following matrix
    fro = np.array([[0, 1], [0, 1], [1, 0]])
    # number row corresponds to the soil emiting compartment
    soilSurf_dic = {
        "Impacted_Soil_Surface": 0,
        "Background_Soil_Surface": 1,
        "Beaches_Soil_Surface": 2,
    }
    # column number corresponds to the recieving compartment

    # In this example of fdd all runoff goes to surface freshwater. To be discussed later

    k_runoff = runoff_rate * fro[soilSurf_dic[particle.Pcompartment.Cname]]
    k_runoff = k_runoff.tolist()

    return k_runoff


def beaching(particle):
    # Transport from surface coastal water to background soil surface.
    # We assume that beaching rate is 1/30 of the transport rate of plastic to open ocean based on https://doi.org/10.1038/s41561-023-01216-0

    if particle.Pcompartment.Cname == "Coast_Surface_Water":

        k_adv = float(particle.Pcompartment.waterFlow_m3_s) / float(
            particle.Pcompartment.Cvolume_m3
        )
        k_beaching = (1 / 30) * k_adv
    else:
        k_beaching = 0

    return k_beaching


def wind_trasport(particle):
    # diffusive transport of particles via wind speed (we should not need this process since ther is onlt one air compartment)
    # to be formulated as funcion of compartment property: wind_speed_m_s
    k_wind_transport = 0
    return k_wind_transport


def dry_deposition(particle, dict_comp):
    # particles depossition from air to soil or water compartments

    # CORRECT NAMING in all the code!!

    # Discuss if to use the dry depossition fractions of distribution here or move it into the fill_interactions function as done for runoff and fragments (we would contruct a dry deposition distribution matrix with the corresponding surface area ratios)

    # Based on figure 6.4 in the Handbook of Chemical Mass Transport in the Environment (2011).

    dd_rate = 7.91e-6

    dd_rate_dict = {
        "e": dd_rate * 1000,
        "d": dd_rate * 100,
        "c": dd_rate,
        "b": dd_rate / 100,
        "a": dd_rate / 1e4,
    }

    k_dry_depossition = [
        dd_rate_dict[particle.Pcode[0]]
        * (
            float(dict_comp[c].CsurfaceArea_m2)
            / float(dict_comp["Air"].CsurfaceArea_m2)
        )
        for c in list(dict_comp.keys())
        if "Surface" in c
    ]
    return k_dry_depossition


def wet_deposition(particle, dict_comp):
    # Currently turned off
    # particles depossition from air to soil or water compartments via rainfall
    # wont be formulated as function of rainfall intensity but dependent on the average rain events per year. we asume that any rain event will trigger the depossition of the particles regardless of rainfall intensity
    # IN  SimpleBox for Plastics rate constant wet depossition 1.17E-1(s-1) Has to be corrected by the number of wet event and duration...so mean rate of depossition will be used
    # wd_rate=?
    # k_dry_depossition = wd_rate*float(particle.Pcompartment.CsurfaceArea_m2)/float(dict_comp["Air"].CsurfaceArea_m2)

    k_wet_depossition = 0
    return k_wet_depossition


def sea_spray_aerosol(particle):
    # particles resuspension from ocean and coastal surface waters to air
    # REF: BETR global approach for ocean-air resuspension (10^-10 m/h)

    ssa_rate = 10e-10 / 60 / 60

    k_sea_spray_aerosol = ssa_rate / float(particle.Pcompartment.Cdepth_m)

    ssa_rate_dict = {
        "a": ssa_rate,
        "b": ssa_rate,
        "c": ssa_rate,
        "d": ssa_rate / 100,
        "e": ssa_rate / 1e4,
    }

    k_sea_spray_aerosol = ssa_rate_dict[particle.Pcode[0]] / float(
        particle.Pcompartment.Cdepth_m
    )

    # k_sea_spray_aerosol = 0

    return k_sea_spray_aerosol


def sequestration_deep_soils(particle):

    # From The OECD tool: MTC3sink = 0.05 * MTCsconv (m/h)soil solids convection to the center of the earth. MTCsconv = 4.54 * 10 ^-7 (m/h)'soil side solid phase convection MTC

    MTCsconv = 4.54e-7

    # K_burial=MTC3sink (m/s) *SA (m2)/V(m3)

    k_sequestration_deep_soils = (0.05 * MTCsconv / (60 * 60)) / float(
        particle.Pcompartment.Cdepth_m
    )

    return k_sequestration_deep_soils
