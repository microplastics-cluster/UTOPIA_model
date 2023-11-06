# UTOPIA_model

A multimedia unit world open-source model for microplastics 

The model comprises 17 compartments and tracks the fate of multiple sizes (5 size bins covering the a size range from nano- to milimeters in size) and four different speciation states of microplastic though solving an overall mass balance that is defined by a system of coupled first-order differential equations.  The masses, m (kg), and or particle number (N) of the different microplastic forms and sizes are obtained as the steady-state solutions of the mass balance equations for all compartments.

![image](https://user-images.githubusercontent.com/58487662/188823913-03dd5f50-2b2c-445d-914a-913637d7bf5b.png)

![image](https://user-images.githubusercontent.com/58487662/188824142-892a10e0-ec4c-42af-adfc-a6a626a35808.png)

The UTOPIA model is being developed based on experiences and knowledge acquired from a previous project ECO48 project Nano2Plast: Extending nanoparticle models to open source models of the fate and transport of plastic in aquatic systems. Therefore, the aquatiac processess included in The Full Multi are used within this project. Further processess of transport between non-aquatic compartments such as sea spray aerosol resuspension to air, runoff of plastics from land, dry and wet depossition of plastics into surface compartments, have been included in UTOPIA (although pending of definition).

![image](https://user-images.githubusercontent.com/58487662/188827636-eeab9b13-f9e1-4e3e-91aa-8041f5b31ba9.png)

1.1. Advective transport 
1.2. Settling 
1.3. Rising 
1.4. Sediment resuspension 
1.5. Diffusion
1.6. Dry deposition 
1.7. Wet deposition/scavenging
1.8. Runoff transport
1.9. Percolation 
1.10. Bioturbation + Tilling
1.11. Sea spray aerosol 
1.12. Soil to air resuspension
1.13. Mixing ??
1.14. Burial in the sediments
1.15. Sequestration within deep soils


## Model Workflow

Description of the workflow followed to develope UTOPIA

### 1- Select input parameters through interactive voila dashboards to generate input csv files 
(to be done)

### 2- Generate model objects by reading on input files:

  -Model boxes (only one box for UTOPIA)
  
  -Model compartments (17 compartments for UTOPIA)
  
  -Particles (one object per MPform (e.g. freeMP, heteroaggregated, biofouled and heteroaggregted and biofouled) and per size fraction (defined 5 size fractions in the range of 0.5 um to 5mm separated by a factor of 10))

### 3-Connect objects to form the modelling system:

 -Assign compartments to boxes
  
 -Add particles to compartments
 
 -Associate particles to the compartments
  
### 4-Parameterise concentrations/emissions

(to be done)

### 5-Generate model processes input parameters table based on the stablished model structure (key parameters such as attachment efficiency, degradation times, fragmentation times, etc.)

### 6-Estimate processes rate constants (per particle in each compartment and model box)

### 7-Generate system of differential equations:

![image](https://user-images.githubusercontent.com/58487662/186609599-c75bb341-45f4-4bf4-a055-fb332aff3756.png)

    -Generate Matrix of interactions (MFullMulti)
    -Build system of differential equations

### 8-Solve mass balance (in steady state for UTOPIA)

### 9-Results presentation through interactive voila dashboards?

  -Mass and particle number concentrations
  -Mass and numbers distribution as fraction of the total mass or total number
  -Exposure metrics (Overall persistence, characteristic travel distance, transfer efficiency)

(to be done)
