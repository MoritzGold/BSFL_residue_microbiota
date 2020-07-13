# BSFL_residue_microbiota
Data and analyses for the manuscript "Identification of high-abundant bacteria in black soldier fly larvae rearing residues"

# waste_composition.R
code for analyses of the substrate nutrient composition

### data = waste_composition.xlsx

**Diet:** household_waste or canteen_waste  
**Replicate:** analyses replicate per diet and parameter  
**Parameter:** acid detergent fibre (ADF), acid detergent fibre (NDF), crude fibre (CF), ash, volatile solids (=organic matter), nitrogen, lipids  
**Unit:** perdm, % DM of freeze-dried sample, already corrected for residual moisture  
**Value:** analysis result

### data = non_fibre_carbohydrates.xlsx

**Sheet1:** fructose (see Available Carbohydrates K-ACHDF, Megazyme, Wicklow, Ireland)

**sample:** household_waste, canteen_waste, blank (=water)  
**blank_ref:** 1 or 2, associated the sample measurment to the respective blank measurement  
**Analysis_replicate:** technical analyses replicate for same sub-sample  
**Replicate:** replicate analyes for different sub-samples  
**sample_weight_g:** sample weight used for the analyis in g  
**final_volume_mL:** final volume in which the absorbance was determined  
**A1:** absorbance 1  
**A2:** absorbance 2  
**A3:** absorbance 3  
**moisture_content_perc:** moisture content of the freeze-dried sample in % DM

**Sheet1:** glucose (see D-Glucose GOPOD K-GLUC, Megazyme, Wicklow, Ireland)

**sample:** household_waste, canteen_waste, blank (=glucose standard) 
**blank_ref:** 1 or 2, associated the sample measurment to the respective control measurement  
**Replicate:** replicate analyes for different sub-samples  
**delta_A_1:** absorbance reading 1 (relative to blank=water)  
**delta_A_2:** absorbance reading 2 (relative to blank=water)  
**volume_mL:** final volume in which the absorbance was determined  
**weight_mg:** sample weight used for the analyis in mg  
**moisture_content_perc:** moisture content of the freeze-dried sample in % DM

**Sheet3:** starch (see Total Starch Assay K-TSTA, Megazyme, Wicklow, Ireland)  
**sample:** household_waste, canteen_waste, blank (=glucose standard)  
**blank_ref:** 1 or 2, associated the sample measurment to the respective control measurement  
**Replicate:** replicate analyes for different sub-samples  
**delta_A_1:** absorbance reading 1 (relative to blank=water)  
**delta_A_2:** absorbance reading 2 (relative to blank=water)  
**weight_mg:** sample weight used for the analyis in mg  
**moisture_content_perc:** moisture content of the freeze-dried sample in % DM
  
# treatment_performance.R
code for analyses of the rearing performance indicators

### data = treatment_performance.xlsx

**Sample_ID:** sample identifier, syntax: xxyzz, xx=Day, y=Diet, zz=Replicate
**Day:** measurement day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment   
**Type:** Sample (with larvae) or Control (without larvae)  
**Diet:** C=canteen waste, S=sterile canteen waste, H=household waste  
**Replicate:** biological replicate per treatment (combination of day, biowaste, type)  
**Number_larvae_in:** number of larvae placed in container at the beginning of the experiment  
**Number_Larvae:** number larvae manually counted after removal of the container at this measurement day   
**DM_Larvae_mg:** dry mass per larvae, in mg per larva    
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

# phyloseq.R
Code for metagenomic analyses with phyloseq and vegan. Uses psychem from physical_chemical_temp.R and performance_messy from treatment_performance.R.

### metedata = code/p529_run190225_MapFile.txt

**SampleID1:** Sample identifier. Syntax xxyzw - x=Day, y=Diet, z=Type, w=Replicate.  
**Diet:**  Diet used for larval rearing. C=canteen waste, H=household waste, S=sterile canteen waste.  
**Day:** Rearing day. 0, 3, 6, 9, 12. 0 marks the beginning of the experiment with larvae having a weight of 0.5 mg DM. Day 12 marks harvesting of larvae from the residue.
**Type:**  Larvae, residue, control or mock samples.  
**Type_detail:**  More detailed description of type.  
**Diet_detail:** Abbreviation of diet.  
**Others:** Full description diet.  
**BSFL_detail:** Denotes whether samples was with ot without larvae.  
**Replicate:** Replicate of diet, day and sample type combination.  
  
# ampvis.R
Code for metagenomic analyses with the ampvis2 package and comparison of our results to the literature

### data = wynants.xlsx
Data from Wynants, E., Frooninckx, L., Crauwels, S., Verreth, C., De Smet, J., Sandrock, C., et al. (2018). Assessing the Microbiota of Black Soldier Fly Larvae (Hermetia illucens) Reared on Organic Waste Streams on Four Different Locations at Laboratory and Large Scale. Microb. Ecol. 

### data = Bruno.xlsx
Data from Bruno, D., Bonelli, M., De Filippis, F., Di Lelio, I., Tettamanti, G., Casartelli, M., et al. (2019). The intestinal microbiota of Hermetia illucens larvae is affected by diet and shows a diverse composition in the different midgut regions. Appl. Environ. Microbiol. 85, e01864-18. 
