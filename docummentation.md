# UTOPIA documentation 
### OBJECTS
•	Model box (UTOPIA)
•	Model compartments
•	Particulates: free, BF and SPM

##### Box Model attributes:
Bname   
Bdepth_m    
Blength_m   
Bwidth_m    
Bvolume_m3  
Bconexions  
Compartments =[]    

•	Model compartments

##### Shared compartment class attributes:
Cname   
Cdepth_m    
Clength_m   
Cwidth_m    
Cvolume_m3  
CsurfaceArea_m2 
particles = {"freeMP": [],"heterMP": [],"biofMP": [],"heterBiofMP":[]}Each key corresponds to another dicttionary of size bins    
connexions = []
processess = [
            "degradation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
        ]

##### Shared compartment class functions
assign_box(self, Box)
dd_particles(self, particle)
calc_volume(self)
calc_particleConcentration_Nm3_initial(self)

Still to be inmplemented: 
assign_particlesEmiss(self, emissionsFile)
assign_backgroundMPConc()


##### Compartment subclasses with specific attributes:
Subclasses (inheritances) of the class compartment add extra attributes to the compatment that define the type of compartment

UTOPIA_surfaceSea_water_compartments = ["Ocean_Surface_Water", "Coast_Surface_Water"]

UTOPIA_water_compartments = [
    "Ocean_Mixed_Water",
    "Ocean_Column_Water",
    "Coast_Column_Water",
    "Surface_Freshwater",
    "Bulk_Freshwater",
]

UTOPIA_deep_soil_compartments = [
    "Beaches_Deep_Soil",
    "Background_Soil",
    "Impacted_Soil",
]

UTOPIA_soil_surface_compartments = [
    "Beaches_Soil_Surface",
    "Background_Soil_Surface",
    "Impacted_Soil_Surface",
]

UTOPIA_sediment_compartment = [
    "Sediment_Freshwater",
    "Sediment_Ocean",
    "Sediment_Coast",
]

UTOPIA_air_compartments = ["Air"]

•	compartment_water specific attributes:
        SPM_mgL,
        waterFlow_m3_s,
        T_K,
        G,
        flowVelocity_m_s
        waterFlow_m3_s

•	compartment_surfaceSea_water specific attributes:
    extra processess: "sea_spray_aerosol","beaching"

•	compartment_sediment specific attributes:
    extra processess = "sediment_resuspension","burial"

•	compartment_soil_surface specific attributes:
    extra processess = "runoff_transport","percolation",
            "soil_air_resuspension","soil_convection".

•	compartment_deep_soil specific attributes:
    extra processess =  "sequestration_deep_soils", "soil_convection"

•	compartment_air specific attributes:
    wind_speed_m_s
    I_rainfall_mm
    extra processess = "wind_trasport","dry_deposition", "wet_deposition"



### MODEL ASUMPTIONS

#### Processess parmeterizations and asumptions

•	Discorporation

Defined as loos of corporeal nature of the particle (not a particle anymore and cosidered as a loss of the particle from the system)

The discorporation rate is assumend to be size, aggregation state and compartment dependent. Values for discorporation rates in the different size fractions, aggregation astates and compartments are estimated from an stabished value of discorporation for free 50 um microplastic particles in surface water (t_half_deg_free), wich is multiplied by a size, aggregation state and compartment factor:

    -Size:  the higuer the surface area to volume ratio of the particle the fastest the process happens according to our assumprions. Therefore, the smaller particles fragment at faster rates. size_factor = (d_base**2) / (d_frac**2), where d_base= 50um
    -Aggregation state asumptions: Heteroaggregated particles degrade 10 times slower than the free MPs (heter_deg_factor = 10) and Biofouled particles degrade 2 times faster than the free MPs (biof_deg_factor= 1/2)
    - Compartment assumptions: we assume that in the surface water compartments both degradation and fragmentation are fastest, in soil surface and deeper water compartments both rates are 10 times slower (factor_deepWater_soilSurface = 10) and in sediments and deeper soil compartments they both are 100 times slower (factor_sediment = 100)

•	Fragmentation

Modelled as a size-dependent process based on an estimated rate constant (𝑘frag_gen= 1/tfrag_gen_d) for fragmentation of pristine (free) particles in the largest (x=5000μm => mp5 => e) size class and in the ater surface compartments.

The fragmentation timescales are deteremined from the stablished fragmentation half time of 36.5 days for the biggest size fraction in free form in the surface water compartments following the parameters chosen in Domercq et al. 2021.

Assumptions: 
Fragmentation is cosidered aggreagtion state and compartment dependent:
    -Aggregation state:fragmentation of the heteroaggregated MPs is stablished as being 100 slower than fragmentation of the Free MPs (heter_frag_factor = 100) and fragmentation of the biofouled particles as being 2 times slower (biof_frag_factor = 2).
    -Compartments: we assume that in the surface water compartments both degradation and fragmentation are fastest, in soil surface and deeper water compartments both rates are 10 times slower (factor_deepWater_soilSurface = 10) and in sediments and deeper soil compartments they both are 100 times slower (factor_sediment = 100)

 The fragmentation relation between size bins is stablished using fragment size distribution matrix (https://microplastics-cluster.github.io/fragment-mnp/advanced-usage/fragment-size-distribution.html). The values of the matrix is determined by the fragmentation inder FI that selects a vlaue between 0 and 1 representing each extream a fragmentation syle as follows:

 - Erosive fragmentation (FI=0): In this scenario the particles are being eroded on their surface and therefore most of their mass remain in their same size fraction and samall fraction in going to the samllest size bins. Its representative fsd is:

         [[0, 0, 0, 0, 0],

         [1, 0, 0, 0, 0],

         [0.99, 0.01, 0, 0, 0],

         [0.999, 0, 0.001, 0, 0],

         [0.9999, 0, 0, 0.0001, 0],]

- Sequential fragmentation (FI=1): in this scenario each size fraction breacks down completely into the next smallest size bin.
Its representative fsd is:

         [[0, 0, 0, 0, 0],

         [1, 0, 0, 0, 0],

         [0, 1, 0, 0, 0],

         [0, 0, 1, 0, 0],

         [0, 0, 0, 1, 0],]

By choosing a value between 0 and 1 the user can select a fragmentation style in between both extremes. (i.e. FI=0.5 will represent the mixed fragmentation style)

•	Heteroaggregation: formulated as in The FullMulti 
•	Heteroaggregate brackup: formulated as in The FullMulti, asumed 10E8 times slower than the formation of heteroaggregates.

•	Advective_transport (same as in The Full Multi)
•	Mixing (references in the RC_generator.py file)
•	Settling (same as in the Full Multi, needs to be updated!)
•	Rising (same as in the Full Multi, needs to be updated!)
•	Biofouling (same as in the Full Multi)
•	Sediment resuspension. Currently placeholder values. To be revisited
•	Burial. Currenlty place holder values. To be revisited
•	Soil air resuspension: We should include a density factor...
•	Soil convection: mixing of soil particles via bioturbation and freeze/thaw cycles
•	Percolation: downwards movement of particles in soil via infiltrated water. To be defined/formulated! (currently =0)
•	Runoff transport: from BETR global approach and mass of runoff distributed to the recieving compartments defined by fro matrix. In this example of fro all runoff goes to surface freshwater. To be reviewed!!
•	Beaching: Transport from surface coastal water to background soil surface. We assume that beaching rate is 1/30 of the transport rate of plastic to open ocean based on https://doi.org/10.1038/s41561-023-01216-0 
•	Dry_deposition: 

## MODEL MODIFICATIONS

•	COMPARTMENTS    

To add/change UTOPIA compartments go to compartments subclasses python file and add new compartment to the compartment list or water, soil or air compartment.
Add the compartment to the imput_compartments.csv file
If new attributes are added to the class this have to also be added in the class definition in the compartmentsSubclassess.py file as None attribute.

•	Rate constants


Time limit can be included in the rate constants using the function timeLimit_particles_RC(system_particle_object_list, k), where k is the maximum value that k can take corresponding to it time limit 1/tlim (currently 30min on the processes that exceeds that speed (k > 0.000556)))


