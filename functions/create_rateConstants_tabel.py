# function to generate rate contants table and save as a .csv file
import pandas as pd

from objects.compartment import Compartment

processList = [
    "k_discorporation",
    "k_fragmentation",
    "k_heteroaggregation",
    "k_heteroaggregate_breackup",
    "k_settling",
    "k_rising",
    "k_advective_transport",
    "k_mixing",
    "k_biofouling",
    "k_sediment_resuspension",
    "k_burial",
    "k_sediment_transport",
    "k_defouling",
]


def create_rateConstants_table(system_particle_object_list, fileName):
    df_dict = {
        "River Section": [],
        "Compartment": [],
        "MP form": [],
        "Size Bin": [],
        "Rate Constants": [],
    }

    for p in system_particle_object_list:
        df_dict["River Section"].append(p.Pcompartment.CBox.Bname)
        df_dict["Compartment"].append(p.Pcompartment.Cname)
        df_dict["MP form"].append(p.Pform)
        df_dict["Size Bin"].append(p.Pname[:3])
        df_dict["Rate Constants"].append(p.RateConstants)

    df = pd.DataFrame(df_dict)
    df2 = df["Rate Constants"].apply(pd.Series)
    df = df.drop(columns="Rate Constants")
    df3 = pd.concat([df, df2], axis=1)
    df4 = df3.fillna(0)

    df4.to_csv(fileName, index=False)

    return df3
