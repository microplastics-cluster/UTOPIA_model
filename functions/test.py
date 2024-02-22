from comp_loop_model_massBalance import model_run

# Model set up

# Size fraction:
# a= 0.5 um
# b= 5 um
# c= 50 um
# d= 500 um
# e= 5000 um
size_dict = dict(zip(["a", "b", "c", "d", "e"], [0.5, 5, 50, 500, 5000]))

# Aggregation state:
# A= Free MP
# B= heteroaggregatedMP
# C= biofouled MP
# D= biofouled and heteroaggregated MP
particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))

MP_form_dict_reverse = {v: k for k, v in particle_forms_coding.items()}

## Compartments:
# 0= 'Ocean_Surface_Water'
# 1= 'Ocean_Mixed_Water'
# 2= 'Ocean_Column_Water'
# 3= 'Coast_Surface_Water'
# 4= 'Coast_Column_Water'
# 5= 'Surface_Freshwater'
# 6= 'Bulk_Freshwater'
# 7= 'Sediment_Freshwater'
# 8= 'Sediment_Ocean'
# 9= 'Sediment_Coast'
# 10= 'Urban_Soil_Surface'
# 11= 'Urban_Soil'
# 12= 'Background_Soil_Surface'
# 13= 'Background_Soil'
# 14= 'Agricultural_Soil_Surface'
# 15= 'Agricultural_Soil'
# 16= 'Air'


model_run(
    mp_imputFile_name="\inputs_microplastics.csv",
    comp_impFile_name="\inputs_compartments.csv",
    comp_interactFile_name="\compartment_interactions.csv",
    q_mass_g_s=1,
    size_bin="e",
    compartment="Ocean_Surface_Water",
    MP_form="freeMP",
)
