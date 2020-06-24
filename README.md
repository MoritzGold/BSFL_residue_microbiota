# BSFL_residue_microbiota
Data and analyses for the manuscript "Identification of high-abundant bacteria in black soldier fly larvae rearing residues"

# treatment_performance.R
code of analyses of the rearing performance indicators

### data = treatment_performance.xlsx

**Sample_ID:** sample identifier, syntax: xxyzz, xx=Day, y=Diet, zz=Replicate
**Day:** measurement day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment   
**Type:** Sample (with larvae) or Control (without larvae)  
**Diet:** C=canteen waste, S=sterile canteen waste, H=household waste  
**Replicate:** biological replicate per treatment (combination of day, biowaste, type)  
**Number_larvae_in:** number of larvae placed in container at the beginning of the experiment  
**Number_Larvae:** number larvae manually counted after removal of the container at this measurement day   
**DM_Larvae_mg:** dry mass per larvae  
**N_Larvae_percDM:** nitrogen content of the freeze-dried larval biomass, in % DM  
**C_Larvae_percDM:** carbon content of the freeze-dried larval biomass, in % DM  
**Ash_Larvae_percDM:** ash content of the freeze-dried larval biomass, in % DM  
**VS_Larvae_percDM:** volatile solids (=organic matter) content of the freeze-dried larval biomass, in % DM  
**Feed_in_gDM:** gram of feed provided, in g DM  
**Feed_out:** residue recovered, in g DM  
**Waste_reduction_DM:** waste reduction, calculated according to Gold, M., Cassar, C. M., Zurbrügg, C., Kreuzer, M., Boulos, S., Diener, S., et al. (2020b). Biowaste treatment with black soldier fly larvae: Increasing performance through the formulation of biowastes based on protein and carbohydrates. Waste Manag. 102, 319–329.  

# counts.R
code for analyses of TVC, LAB and fungi counts in the substrate, residue and larvae

### data = counts.xlsx

**day**: measurement day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment.  
**biowaste**: C=canteen waste, S=sterile canteen waste, H=household waste  
**type**:larva, residue or no_larvae (control without larvae)  
**replicate**: biological replicate per treatment (combination of day, biowaste, type)  
**parameter**: describes the microbial number parameter, either TVC, LAB (=lab) or fungi (=fy) counts  
**value:** TVC, LAB and fungi count  

# physical_chemical_temp.R
code for analyses of physio-chemical parameters and temperature in the residue  
uses performance_psychem from performance_psychem.R

### data = physico_chemical.xlsx

**day:** measurement day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment     
**biowaste**: C=canteen waste, S=sterile canteen waste, H=household waste   
**type:** residue with larvae or without larvae  
**replicate:**  replicate per treatment (combination of day, biowaste, type)  
**parameter:** ash (100%-VS %DM), aw (water activity), C (nitrogen content %DM), MC (moisture content %DM), N (nitrogen content %DM), pH, VS (volatile solids = organic matter %DM)  
**analyses_replicate:** technical analyses replicate (repeated measurment of the parameter for the same replicate)  
**value** measurment results of the parameter  

### data = temp.xlsx

**day:** Temperature recordings between day 0-3 (3), 3-6 (6), 6-9 (9) or 9-12 (12)  
**biowaste:** type of waste (S=sterile canteen waste, C = canteeen waste, H = household waste), CC = temperature recording
in the climate chamber.  
**type:** larvae or no larvae (controls without larvae)  
**replicate:** replicate per treatment (combination of day, biowaste, type)  
**parameter:** temp. (=temperature)  
**value:** temperature in °C  
