import os
import pandas as pd

# print("getcwd: ", os.getcwd())
# print("__file__: ", __file__)
# print("basename: ", os.path.basename(__file__))
# print("dirname: ", os.path.dirname(__file__))

target_path = os.path.join(
    os.path.dirname(__file__), "../inputs/processInputs_table.csv"
)

print(target_path)

process_inputs_df = pd.read_csv(filepath_or_buffer=target_path)

print(process_inputs_df)
