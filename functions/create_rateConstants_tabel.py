# function to generate rate contants table and save as a .csv file
import pandas as pd

from objects.compartment import Compartment


def create_rateConstants_table(system_particle_object_list):
    df_dict = {
        "Compartment": [],
        "MP_form": [],
        "Size_Bin": [],
        "Rate_Constants": [],
    }

    for p in system_particle_object_list:
        df_dict["Compartment"].append(p.Pcompartment.Cname)
        df_dict["MP_form"].append(p.Pform)
        df_dict["Size_Bin"].append(p.Pname[:3])
        df_dict["Rate_Constants"].append(p.RateConstants)

    df = pd.DataFrame(df_dict)
    df2 = df["Rate_Constants"].apply(pd.Series)
    df = df.drop(columns="Rate_Constants")
    df3 = pd.concat([df, df2], axis=1)

    return df3
