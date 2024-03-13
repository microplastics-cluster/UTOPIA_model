import pandas as pd
import os


def write_MPinputs_table(
    MPdensity_kg_m3,
    MP_composition,
    shape,
    N_sizeBins,
    big_bin_diameter_um,
    runName,
    inputs_path,
):
    size_distribution = [big_bin_diameter_um]
    N = N_sizeBins - 1
    for i in range(N):
        Bsize = size_distribution[-1]
        newSize = Bsize / 10
        size_distribution.append(newSize)
    size_distribution.reverse()
    if shape == "sphere":
        data = {
            "Name": ["mp" + str(i + 1) for i in range(N_sizeBins)],
            "form": ["freeMP"] * N_sizeBins,
            "composition": [MP_composition] * N_sizeBins,
            "shape": [shape] * N_sizeBins,
            "density_kg_m3": [MPdensity_kg_m3] * N_sizeBins,
            "dimensionX_um": [i / 2 for i in size_distribution],
            "dimensionY_um": [i / 2 for i in size_distribution],
            "dimensionZ_um": [i / 2 for i in size_distribution],
        }

    else:
        print("Shape not supported yet")
    df = pd.DataFrame(data)

    inputsMP_fileName = os.path.join(
        inputs_path, "inputs_microplastics" + runName + ".csv"
    )

    df.to_csv(inputsMP_fileName, index=False)

    return inputsMP_fileName
