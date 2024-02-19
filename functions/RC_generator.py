# Extension of RS_generator module from the FUll Multi containing functions to calculate all rate constants


###So far programmed so that if a parameter value for estimating a rate constant is missing (i.e. NOne or nan), the rate constant is equal to zero....this has to be fixed!! Example,when no depth ogf a compartment is given the rate constant of burial or resuspension is equal to 0

import math
import pandas as pd
import os
import numpy as np


# import file storing required constants
from helpers.globalConstants import *

# Read input data file

process_inputs_df = pd.read_csv(
    filepath_or_buffer=os.path.join(
        os.path.dirname(__file__), "../inputs/processInputs_table.csv"
    )
)


def discorporation(particle):
    # degradation estimations
    # discorporation state will be given as output after runing
    # the model with no discorporation. degradation state will be given in time units as residence time in the compartment

    # Change process name from degradation
    # to discorporation from corporeal

    """relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    # degradation half-life of MPs used as input is in days
    cond = (
        (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
        & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
    )
    t_half_d = float(process_inputs_df.loc[cond, "thalf_deg_d"])

    # degradation rate constant
    k_deg = math.log(2) / (t_half_d * 24 * 60 * 60)

    return k_deg


def fragmentation(particle):

    # modelled as a size-dependent process based on an estimated rate constant (𝑘frag_gen= 1/tfrag_gen_d)
    # for fragmentation of pristine particles in the largest (x=5000μm => mp1 => e) size class.

    # estimate fragmentation relation between size bins using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html)

    # Fragmentation of heteroaggregated particles is assumed negligible in the default model formulation

    cond = (
        (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
        & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == "mp1")
    )
    t_frag_d = process_inputs_df.loc[cond, "tfrag_gen_d"].item()

    if t_frag_d == "NAN":
        # It is assumed that fragmentation doesnt occur in heteroaggregated particles (heterMP or heterBiofMP)
        frag_rate = 0
        fragments_formed = 0
    else:
        if (
            particle.Pname[0:3] == "mp5"
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
    # each article fractions into fragments of samller sizes and the distribution is expresses via the fragment size distribution matrix fsd. # In this matrix the smallest size fraction is in the first possition and we consider no fragmentation for this size class
    size_dict = {chr(i): i - ord("a") for i in range(ord("a"), ord("e") + 1)}
    # fsd = np.array(
    #     [
    #         [0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 0],
    #     ]
    # )

    fsd = np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0.5, 0.5, 0, 0, 0],
            [0.6, 0.2, 0.2, 0, 0],
            [0.7, 0.15, 0.1, 0.05, 0],
        ]
    )
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


def heteroaggregation(particle, spm):
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
            (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
            & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
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


# k_hetAgg


def heteroaggregate_breackup(particle, spm):
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
                (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
                & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
                & (process_inputs_df["MPform"] == "freeMP")
                & (process_inputs_df["sizeBin"] == particle.Pname[0:3])
            )
        elif particle.Pform == "heterBiofMP":
            cond_alpha = (
                (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
                & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
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


# k_aggBreakup


def advective_transport(particle):
    k_adv = float(particle.Pcompartment.waterFlow_m3_s) / float(
        particle.Pcompartment.Cvolume_m3
    )

    # if particle.Pcompartment.waterFlow_m3_s != "nan":
    #     k_adv = float(particle.Pcompartment.waterFlow_m3_s) / float(particle.Pcompartment.Cvolume_m3)
    # else:
    #     k_adv = float(particle.Pcompartment.flowVelocity_m_s) * (
    #         float(particle.Pcompartment.Cdepth_m)
    #         * float(particle.Pcompartment.Cwidth_m)
    #         / float(particle.Pcompartment.Cvolume_m3)
    #     )
    # advective transport

    # Based on Praetorius et al. 2012: Kflow = v_riv_flow*(Aw1/Vw1)
    # Being v_riv_flow the river flow velocity in ms-1, Aw1 is the crossectional
    # area of the flowing water and Vw1 the volume of the box of moving water.

    """So far advection was described for a river with changing discharges, TO BE ADAPTED FOR THE UTOPIA"""

    return k_adv


def mixing(particle, dict_comp):
    # Now adapted to UTOPIA's compartments and changed rates-- In progress--
    # k_mix has to be multiplied by the compartment volume ratio calculated with the interacting compartment volume

    k_mix_up = (
        10**-2
    )  # (1): <Handbook of Chemical Mass Transport in the Environment> Edited by Louis J. Thibodeaux, Donald Mackay (DOI: 10.1201/b10262)

    k_mix_down = (
        10**-3
    )  # (2): <Handbook on Mixing in Rivers> Edited by J.C. Rutherford (Water and Soil Miscellaneous Publication No. 26. 1981. 60pp.ISSN 0110-4705)

    if particle.Pcompartment.Cname == "Ocean_Mixed_Water":
        k_mix = [
            k_mix_up
            * float(dict_comp["Ocean_Surface_Water"].Cvolume_m3)
            / float(particle.Pcompartment.Cvolume_m3),
            k_mix_down
            * float(dict_comp["Ocean_Column_Water"].Cvolume_m3)
            / float(particle.Pcompartment.Cvolume_m3),
        ]
        # {"mix_up": k_mix_up, "mix_down": k_mix_down}

    elif particle.Pcompartment.Cname in [
        "Ocean_Column_Water",
        "Coast_Column_Water",
        "Bulk_Freshwater",
    ]:
        connecting_comp = [
            key
            for key, value in particle.Pcompartment.connexions.items()
            if (isinstance(value, list) and "mixing" in value)
        ][0]

        k_mix = (
            k_mix_up
            * float(dict_comp[connecting_comp].Cvolume_m3)
            / float(particle.Pcompartment.Cvolume_m3)
        )
    elif particle.Pcompartment.Cname in [
        "Ocean_Surface_Water",
        "Coast_Surface_Water",
        "Surface_Freshwater",
    ]:
        connecting_comp = [
            key
            for key, value in particle.Pcompartment.connexions.items()
            if (isinstance(value, list) and "mixing" in value)
        ][0]
        k_mix = (
            k_mix_down
            * float(dict_comp[connecting_comp].Cvolume_m3)
            / float(particle.Pcompartment.Cvolume_m3)
        )
    else:
        k_mix = 0

    return k_mix


"""To be described for UTOPIA"""


def biofouling(particle):
    cond_biof = (
        (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
        & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
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


# k_biof


def defouling(particle):
    # Defouling = degradation of Biofilm.

    cond_defoul = (
        (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
        & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
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


# for the sediment compartment rate constants for resuspension and
# burial in deep sediment are calculated & degradation rate assigned


def sediment_resuspension(particle):
    # When no depth parameter available assign transfer sediment to water rate taken from SimpleBox for Plastics model
    resusp_dict = {
        "Sediment_Freshwater": 1e-9,
        "Sediment_Coast": 1e-10,
        "Sediment_Ocean": 1e-11,
    }

    k_resusp = resusp_dict[particle.Pcompartment.Cname]

    return k_resusp


def burial(particle):
    # When no depth parameter available assign burail rate taken from SimpleBox for Plastics model
    burial_dict = {
        "Sediment_Freshwater": 2.7e-10,
        "Sediment_Coast": 2.7e-11,
        "Sediment_Ocean": 2.7e-12,
    }

    k_burial = burial_dict[particle.Pcompartment.Cname]

    return k_burial


def soil_air_resuspension(particle):
    # To be formulated
    # default value talen from SimpleBox for Plastics as trasnfer rate soil-air in s-1
    # If it would be compartment dependent we can sue a dictionary with values: now we use the same values for all soil surface compartments
    sa_resusp_dic = {
        "Urban_Soil_Surface": 4.68e-24,
        "Background_Soil_Surface": 4.68e-24,
        "Agricultural_Soil_Surface": 4.68e-24,
    }

    k_sa_reusp = sa_resusp_dic[particle.Pcompartment.Cname]

    return k_sa_reusp


def tillage(particle):
    # vertical mixing of particles in soil via physical mixing of top soil
    # defines transport between soil surface and deeper layers
    # to be defined/formulated
    pass


def percolation(particle):
    # downwards movement of particles in soil via infiltrated water
    # to be formulated

    # k_percol = particle.Pcompartment.infiltration_capacity*particle.Pcompartment.precipitation_rate*(float(particle.Pcompartment.Cvolume_m3)/float(particle.Pcompartment.Cdepth_m))/float(particle.Pcompartment.soilPore_waterVolume_m3)

    k_percol = 0

    return k_percol


def runoff_transport(particle):
    # transport from top soil layers to surface waters ["Coast_Surface_Water","Surface_Freshwater"] via runoff water
    # to be formulated
    runooff_dict = {
        "Urban_Soil_Surface": 4.69e-9,
        "Background_Soil_Surface": 4.69e-9,
        "Agricultural_Soil_Surface": 4.69e-9,
    }
    runoff_rate = runooff_dict[particle.Pcompartment.Cname]

    # The total amount of runoff will be distributed into the recieving compartments according to the following matrix
    fro = np.array([[0, 1], [0, 1], [0, 1]])
    # number row corresponds to the soil emiting compartment
    soilSurf_dic = {
        "Urban_Soil_Surface": 0,
        "Background_Soil_Surface": 1,
        "Agricultural_Soil_Surface": 2,
    }
    # column number corresponds to the recieving compartment

    # In this example of fdd all runoff goes to surface freshwater. To be discussed later

    k_runoff = runoff_rate * fro[soilSurf_dic[particle.Pcompartment.Cname]]
    k_runoff = k_runoff.tolist()

    return k_runoff


def wind_trasport(particle):
    # diffusive transport of particles via wind speed (we should not need this process since ther is onlt one air compartment)
    # to be formulated as funcion of compartment property: wind_speed_m_s
    k_wind_transport = 0
    return k_wind_transport


def dry_depossition(particle, dict_comp):
    # particles depossition from air to soil or water compartments
    # to be formulated
    # Default value taken from SimpleBox for Plastics rate constant dry depossition 2.16E-6 (s-1).
    # Will be reformulated to be made size class and recieving compartment dependent

    # Discuss if to use the dry depossition fractions of distribution here or move it into the fill_interactions function as done for runoff and fragments (we would contruct a dry deposition distribution matrix with the corresponding surface area ratios)

    dd_rate = 7.91e-6
    k_dry_depossition = [
        dd_rate
        * (
            float(dict_comp[c].CsurfaceArea_m2)
            / float(dict_comp["Air"].CsurfaceArea_m2)
        )
        for c in list(dict_comp.keys())
        if "Surface" in c
    ]
    return k_dry_depossition


def wet_depossition(particle, dict_comp):
    # particles depossition from air to soil or water compartments via rainfall
    # wont be formulated as function of rainfall intensity but dependent on the average rain events per year. we asume that any rain event will trigger the depossition of the particles regardless of rainfall intensity
    # Default value taken from SimpleBox for Plastics rate constant wet depossition 1.17E-1(s-1) Has to be corrected by the number of wet event nd duration...so mean rate of depossition will be used
    # wd_rate=?
    # k_dry_depossition = wd_rate*float(particle.Pcompartment.CsurfaceArea_m2)/float(dict_comp["Air"].CsurfaceArea_m2)

    k_wet_depossition = 0
    return k_wet_depossition


def sea_spray_aerosol(particle):
    # paticles resuspension from ocean and coastal surface waters to air
    # Default value taken from SimpleBox for Plastics transfer rate water-air (s-1)
    # to be formulated

    k_sea_spray_aerosol = 2.36e-25
    return k_sea_spray_aerosol


def sequestration_deep_soils(particle):
    # to be formulated
    # Default value taken from SimpleBox for Plastics Removal rate from soil in s-1
    k_sequestration_deep_soils = 2.71e-9

    return k_sequestration_deep_soils
