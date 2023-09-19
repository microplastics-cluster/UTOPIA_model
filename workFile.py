import pandas as pd

import numpy as np


def remove_diag(x):
    x_no_diag = np.ndarray.flatten(x)
    x_no_diag = np.delete(x_no_diag, range(0, len(x_no_diag), len(x) + 1), 0)
    x_no_diag = x_no_diag.reshape(len(x), len(x) - 1)
    return x_no_diag


matrix_no_diag = remove_diag(matrix)

plt.hist(matrix_no_diag)


log_data = np.log10(matrix_no_zero)
plt.hist(log_data)


#########################
# Fuction to create compartment objects from csv file
import os
import csv
from helpers.helpers import *


class Record:
    """Hold a record of data."""


inputs_path = os.path.join(os.path.dirname(__file__), "inputs")

compFile = "\inputs_microplastics.csv"
with open(inputs_path + compFile, "r") as f:
    dictReader_obj = csv.DictReader(f)
    particles = list(dictReader_obj)

    class Record:
        """Hold a record of data."""

    # for item in dictReader_obj:
    #     for field, value in item.items():
    #         setattr(particle_record, field,value)
    #     print(dict(item))
    particlesObj_list = []
    for item in particles:
        particle_record = Record()
        for field, value in item.items():
            setattr(particle_record, field, value)
        particlesObj_list.append(particle_record)


##########
# Create connexions attributes as dictionaries for the different #compartments from the compartmentsInteractions file
def set_interactions(compartments, connexions_path_file="compartmentsInteractions.csv"):
    comp_connex_df = pd.read_csv(connexions_path_file)

    for c in compartments:
        df_c = comp_connex_df[["Compartments", c.Cname]].dropna()
        c.connexions = dict(zip(df_c["Compartments"], df_c[c.Cname]))


dict([(key, value) for key, value in data])
