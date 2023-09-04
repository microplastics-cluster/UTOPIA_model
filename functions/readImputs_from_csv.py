# Instantiate class from csv file: each line in the csv will generate one particle object
from objects.particulates import Particulates
from objects.box import Box
from objects.compartmetSubclasess import *
import csv
import pandas as pd
import numpy as np


def instantiateParticles_from_csv(compFile):
    with open(compFile, "r") as f:
        reader = csv.DictReader(f)
        particles = list(reader)

    particlesObj_list = []
    for p in particles:
        particlesObj_list.append(
            Particulates(
                Pname=p.get("Name"),
                Pform=p.get("form"),
                Pcomposition=p.get("composition"),
                Pdensity_kg_m3=float(p.get("density_kg_m3")),
                Pshape=p.get("shape"),
                PdimensionX_um=float(p.get("dimensionX_um")),
                PdimensionY_um=float(p.get("dimensionY_um")),
                PdimensionZ_um=float(p.get("dimensionZ_um")),
            )
        )

    print(
        f"The free MP particles {[p.Pname for p in particlesObj_list]} have been generated"
    )

    return particlesObj_list


def instantiateBoxes_from_csv(boxFile):
    with open(boxFile, "r") as f:
        reader = csv.DictReader(f)
        boxes = list(reader)

    boxesObject_list = []
    for b in boxes:
        boxesObject_list.append(
            Box(
                Bname=b.get("name"),
                Bdepth_m=float(b.get("depth_m")),
                Blength_m=float(b.get("length_m")),
                Bwidth_m=float(b.get("width_m")),
                Bvolume_m3=b.get("Bvolume_m3"),
                Bconexions=b.get("conexions"),
            )
        )
    return boxesObject_list


def parameteriseRiverSections_from_csv(temp_RS_properties, riverSections):
    RSproperties = pd.read_csv(temp_RS_properties)
    for riverSect in riverSections:
        riverSect.T_K = RSproperties.loc[
            RSproperties["Bname"] == riverSect.Bname, "T_K"
        ].item()
        riverSect.spm_mgL = RSproperties.loc[
            RSproperties["Bname"] == riverSect.Bname, "conc_SPM_mg_L"
        ].item()
        riverSect.Ca_mg_L = RSproperties.loc[
            RSproperties["Bname"] == riverSect.Bname, "Ca_mg_L"
        ].item()
        riverSect.DOC_mg_L = RSproperties.loc[
            RSproperties["Bname"] == riverSect.Bname, "DOC_mg_L"
        ].item()
        riverSect.conexions = RSproperties.loc[
            RSproperties["Bname"] == riverSect.Bname, "DOC_mg_L"
        ].item()


def read_connexions_inputs(inputs_path, input_file, compartmentNames_list):
    df = pd.read_csv(inputs_path + input_file)
    dicts = []
    for i in range(len(compartmentNames_list)):
        df_sub = df[df.Compartment == compartmentNames_list[i]]
        dicts.append(
            {
                compartmentNames_list[i]: [
                    {df_sub.Connexion[x]: df_sub.Process[x]} for x in df_sub.index
                ]
            }
        )
    connexions_dict = {k: v for d in dicts for k, v in d.items()}
    return connexions_dict


def instantiate_compartments(inputs_path_file):
    with open(inputs_path_file, "r") as f:
        reader = csv.DictReader(f)
        compartments = list(reader)

    waterComp_objects = []
    soilComp_objects = []
    airComp_objects = []
    for c in compartments:
        if c["Cname"] in UTOPIA_water_compartments:
            waterComp_objects.append(
                compartment_water(
                    Cname=c.get("Cname"),
                    SPM_mgL=c.get("SPM_mgL"),
                    flowVelocity_m_s=c.get("flowVelocity_m_s"),
                    waterFlow_m3_s=c.get("waterFlow_m3_s"),
                    T_K=c.get("T_K"),
                    G=c.get("G"),
                    Cdepth_m=c.get("Cdepth_m"),
                    Cvolume_m3=c.get("Cvolume_m3"),
                )
            )
        elif c["Cname"] in UTOPIA_soil_compartments:
            soilComp_objects.append(compartment_soil(Cname=c.get("Cname")))

        elif c["Cname"] in UTOPIA_air_compartments:
            airComp_objects.append(compartment_air(Cname=c.get("Cname")))

        else:
            pass

    Comp_objects = waterComp_objects + soilComp_objects + airComp_objects

    print(f"The compartments {[c.Cname for c in Comp_objects]} have been generated")

    return Comp_objects


def instantiate_compartments_from_csv(inputs_path_file):

    # Read csv as dictionarys for each line
    compartments = {}
    with open(inputs_path_file, "r") as f:
        reader = csv.DictReader(f)
        for item in reader:
            compartments[item["Cname"].replace(" ", "")] = dict(item)
    # transform nested dictionaries to objects with attributes
    comp_objects_list = []
    for c in compartments:
        obj = Struct(**compartments[c])
        comp_objects_list.append(obj)

    return comp_objects_list


# define a class to convert a nested dictionary to an Object
class Struct:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            if isinstance(value, dict):
                self.__dict__[key] = Struct(**value)
            else:
                self.__dict__[key] = value


# Create connexions attributes as dictionaries for the different #compartments from the compartmentsInteractions file
def set_interactions(compartments, connexions_path_file):

    with open(connexions_path_file, "r") as infile:
        reader = csv.reader(infile)
        array = []
        for row in reader:
            r = []
            for ele in row:
                if "," in ele:
                    r.append(ele.split(","))
                else:
                    r.append(ele)
            array.append(r)
        comp_connex_df = pd.DataFrame(array)
        comp_connex_df.columns = comp_connex_df[0]
        # comp_connex_df = comp_connex_df.set_index("Compartments")
        comp_connex_df = comp_connex_df.drop(index=[0])
        comp_connex_df.replace("", np.nan, inplace=True)

    # comp_connex_df = pd.read_csv(connexions_path_file)

    for c in compartments:
        df_comp = comp_connex_df[["Compartments", c.Cname]].dropna()
        c.connexions = dict(zip(df_comp["Compartments"], df_comp[c.Cname]))

    print("Connexions have been added")
