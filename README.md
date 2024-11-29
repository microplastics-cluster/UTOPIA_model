# UTOPIA_model

A multimedia unit world open-source model for microplastics 

The model comprises 17 compartments and tracks the fate of multiple sizes (5 size bins covering the a size range from nano- to milimeters in size) and four different speciation states of microplastic though solving an overall mass balance that is defined by a system of coupled first-order differential equations.  The masses, m (kg), and or particle number (N) of the different microplastic forms and sizes are obtained as the steady-state solutions of the mass balance equations for all compartments.

![image](https://github.com/microplastics-cluster/UTOPIA_model/blob/main/UTOPIA_building_blocks.png)

![image](https://user-images.githubusercontent.com/58487662/188824142-892a10e0-ec4c-42af-adfc-a6a626a35808.png)

The UTOPIA model is being developed based on experiences and knowledge acquired from a previous project ECO48 project Nano2Plast: Extending nanoparticle models to open source models of the fate and transport of plastic in aquatic systems. Therefore, the aquatiac processess included in The Full Multi are used within this project. Further processess of transport between non-aquatic compartments such as sea spray aerosol resuspension to air, runoff of plastics from land, dry and wet depossition of plastics into surface compartments, have been included in UTOPIA (although pending of definition).

![image](https://github.com/microplastics-cluster/UTOPIA_model/blob/main/UTOPIA_processes.png)

The processess marked in red are not yet included in UTOPIA.

## Model Workflow

Description of the workflow followed to develope UTOPIA

### 1- Select input parameters
Go to the file script_UTOPIA_user.py and follow the instructions given in the code. Currently it is possible to modify:

  Microplastics properties:

  - MPdensity_kg_m3 
  - MP_composition  (Has to match the defined density)

  Compartment properties:

  To do so one should copy the inputs_compartments.csv file and modify its values without changing the format of the file. Once a new file is generated this can be saved with a new name and the new name should be provided for the comp_impFile_name variable of the script_UTOPIA_user.py

  Plastic weathering properties:

  - Select a fragmentation style by choosing a FI value that goes from 0 to 1 to select a scenario between Erosive (FI0=) and sequential (FI=1) fragmentation as described in the code.
  - Type the fragmentation timescale (given in number of days) in: t_frag_gen_FreeSurfaceWater (default values= 36.5)
  - Type the disintegration timescale (given in number of days) in: t_half_deg_free (default value=66000)
  - Fractors that infuence fragmentation and disintegration according to the environmental compartment and aggregation state are defined with default values as listed:
  - heter_deg_factor = 10
  - biof_deg_factor = 1 / 2
  - factor_deepWater_soilSurface = 10
  - factor_sediment = 100
  - biof_frag_factor = 2
  - heter_frag_factor = 100

  Emission scenario:
  describe the porperties of the emitted plastic particles (particle size and form) and the recieving compartment/s and its flow of emission/s

  - size_bin: choose a size fraction from the size_dict dictionary (if chosen the default settings: a= 0.5 um, b= 5 um, c= 50 um, d= 500 um, e= 5000 um)
  - MP_form: Choose from MPforms_list (freeMP,heterMP,biofMP and heterBiofMP)
  - If emission are targeted to a single compartment the user should define:
    - input_flow_g_s
    - emiss_comp
  - If there are emissions into several compartments the user should add the corresponding input flows per compartment in the dictionary q_mass_g_s_dict (L252). Note that if this second option is chossen the user should double check that the inputs above match the enission scenario targeted.


### 2- Generate model objects by reading on input files:

  -Model boxes (only one box for UTOPIA)
  
  -Model compartments (17 compartments for UTOPIA)
  
  -Particles (one object per MPform (e.g. freeMP, heteroaggregated, biofouled and heteroaggregted and biofouled) and per size fraction (defined 5 size fractions in the range of 0.5 um to 5mm separated by a factor of 10))

### 3-Connect objects to form the modelling system:

 -Assign compartments to boxes
  
 -Add particles to compartments
 
 -Associate particles to the compartments
  
### 4-Parameterise concentrations/emissions


### 5-Generate model processes input parameters table based on the stablished model structure (key parameters such as attachment efficiency, degradation times, fragmentation times, etc.)

### 6-Estimate processes rate constants (per particle in each compartment and model box)

### 7-Generate system of differential equations:

![image](https://user-images.githubusercontent.com/58487662/186609599-c75bb341-45f4-4bf4-a055-fb332aff3756.png)

    -Generate Matrix of interactions (MFullMulti)
    -Build system of differential equations

### 8-Solve mass balance (in steady state for UTOPIA)

### 9-Plot and save results

  -Mass and particle number concentrations
  -Mass and numbers distribution as fraction of the total mass or total number
  -Exposure metrics (Overall persistence, characteristic travel distance, transfer efficiency)
  -Emission fractions


## Instalation guidelines

### Getting started for Windows

Download the repository to your computer by clicking on the green CODE button on the rigth of the repository screen.

### Create, activate, and download dependencies with a virtual environment using venv

### Create a virtual environment named 'venv'
python -m venv venv

### Activate the virtual environment on Windows
venv\Scripts\activate

### Install requirements
pip install -r requirements.txt

### Run server 
```bash
python script_UTOPIA_user.py
```
