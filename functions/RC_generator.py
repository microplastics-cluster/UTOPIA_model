# Extension of RS_generator module from the FUll Multi containing functions to calculate all rate constants

import math
import pandas as pd
import os


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
    # estimate fragmentation relation between size bins (all except smallest size bin)
    # modelled as a size-dependent process based on an estimated rate constant (𝑘frag_gen= 1/tfrag_gen_d)
    # for fragmentation of pristine particles in the largest (1000 μm=mp5) size class.

    # Fragmentation of heteroaggregated particles is assumed negligible in the default model formulation

    cond = (
        (process_inputs_df["modelBox"] == particle.Pcompartment.CBox.Bname)
        & (process_inputs_df["Compartment"] == particle.Pcompartment.Cname)
        & (process_inputs_df["MPform"] == particle.Pform)
        & (process_inputs_df["sizeBin"] == "mp5")
    )
    t_frag_d = process_inputs_df.loc[cond, "tfrag_gen_d"].item()

    if t_frag_d == "NAN":
        # It is assumed that fragmentation doesnt occur in heteroaggregated particles (heterMP or heterBiofMP)
        k_frag = 0
        fragments_formed = 0
    else:
        if (
            particle.Pname[2:3] == "mp1"
        ):  # print("Smallest sizeBin, fragments formed will be considered losses")
            k_frag = k_frag = (
                (1 / (float(t_frag_d) * 24 * 60 * 60))
                * float(particle.diameter_um)
                / 1000
            )
            fragments_formed = 0
        else:
            volume_fragment = (
                4 / 3 * math.pi * (float(particle.radius_m) / 10) ** 3
            )  #!!!only works for bins 10 times smaller!!!
            fragments_formed = float(particle.Pvolume_m3) / volume_fragment
            k_frag = (
                (1 / (float(t_frag_d) * 24 * 60 * 60))
                * float(particle.diameter_um)
                / 1000
            )

    return k_frag  # I have removed the fragments formed from the output to have an homogeneus solution in the table of rate constants (consider dumping this values in another way later when needed (for Mass Balance?))

    # NOTE: to be modified by ECO59


def settling(particle):
    # settling calculations
    """settling can be calculated using different equations (e.g. Stokes,
    modified versions of it or others) or be taken from experimental studies
    !! currently only classical Stokes is implemented (which is probably not
    very realistic and will be updated soon !!"""

    settlingMethod = "Stokes"

    # Settling occurs in all aquatic compartments which should be specified in the comprtment class
    # if particle.Pcompartment.Cname in ["Sediment", "Agricultural Soil","Urban Soil"...]
    #     k_set = 0

    if settlingMethod == "Stokes":
        vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - density_w_21C_kg_m3)
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
        "Ocean Mixed Water",
        "Ocean Column Water",
        "Coast Column Water",
        "Bulk FreshWater",
    ]:
        if settlingMethod == "Stokes":
            vSet_m_s = (
                2
                / 9
                * (float(particle.Pdensity_kg_m3) - density_w_21C_kg_m3)
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

        MP_vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - density_w_21C_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (float(particle.radius_m)) ** 2
        )

        SPM_vSet_m_s = (
            2
            / 9
            * (spm.Pdensity_kg_m3 - density_w_21C_kg_m3)
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


def heteroaggregate_breackup(particle, spm):

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

        MP_vSet_m_s = (
            2
            / 9
            * (float(particle.Pdensity_kg_m3) - density_w_21C_kg_m3)
            / mu_w_21C_kg_ms
            * g_m_s2
            * (float(particle.radius_m)) ** 2
        )
        SPM_vSet_m_s = (
            2
            / 9
            * (spm.Pdensity_kg_m3 - density_w_21C_kg_m3)
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

            k_aggBreakup = 1 / 10 * k_hetAgg
    else:
        k_aggBreakup = 0

    return k_aggBreakup


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
    # Taken from the Full Multi Model, should be adapted to UTOPIA

    if particle.Pcompartment.Cname == "Flowing Water":
        k_mix_up = 10**-10
        k_mix_down = 10**-13
        k_mix = (k_mix_up, k_mix_down)
    elif particle.Pcompartment.Cname == "Surface Water":
        k_mix = (10**-10) * (
            dict_comp["Flowing Water"].Cvolume_m3
            / float(particle.Pcompartment.Cvolume_m3)
        )
    elif particle.Pcompartment.Cname == "Stagnant Water":
        k_mix = (10**-13) * (
            dict_comp["Flowing Water"].Cvolume_m3
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

    k_resusp = 2.3 * 10**-7 / float(particle.Pcompartment.Cdepth_m)

    return k_resusp


def burial(particle):

    k_burial = 5.6 * 10**-7 / float(particle.Pcompartment.Cdepth_m)

    return k_burial


def sediment_transport(particle):
    m_sed_kg = (
        (1 - sed_porosity)
        * sed_density
        * 10**3
        * float(particle.Pcompartment.Cvolume_m3)
    )
    k_sed_trans = v_sed_trans / m_sed_kg

    return k_sed_trans


def soil_air_resuspension(particle):

    k_sa_reusp = 0  # To be formulated

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
    # transport from top soil layers to surface waters via runoff water
    # to be formulated
    k_runoff = 0
    return k_runoff


def wind_trasport(particle):
    # diffusive transport of particles via wind speed
    # to be formulated as funcion of compartment property: wind_speed_m_s
    k_wind_transport = 0
    return k_wind_transport


def dry_depossition(particle):
    # particles depossition from air to soil or water compartments
    # to be formulated
    k_dry_depossition = 0
    return k_dry_depossition


def wet_depossition(particle):
    # particles depossition from air to soil or water compartments via rainfall
    # to be formulated as function of rainfall intensity??
    k_wet_depossition = 0
    return k_wet_depossition


def sea_spray_aerosol(particle):
    # paticles resuspension from ocean and coastal surface waters to air

    # to be formulated

    k_sea_spray_aerosol = 0
    return k_sea_spray_aerosol


def sequestration_deep_soils(particle):
    # to be formulated
    k_sequestration_deep_soils = 0

    return k_sequestration_deep_soils
