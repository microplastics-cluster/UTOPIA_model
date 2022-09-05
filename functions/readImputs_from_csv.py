# Instantiate class from csv file: each line in the csv will generate one particle object
from objects.particulates import Particulates
from objects.box import Box
import csv
import pandas as pd


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
