# BSFL_residue_microbiota
Data and analyses for the manuscript "Identification of high-abundant bacteria in black soldier fly larvae rearing residues"

# counts.R
code on TVC, LAB and fungi counts in the substrate, residue and larvae

### data = counts.xlsx

**day**: measurment day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment.  
**biowaste**: C=canteen waste, S=sterile canteen waste, H=household waste  
**type**:larva, residue or no_larvae (control without larvae)  
**replicate**: replicate per treatment (combination of day, biowaste, type)  
**parameter**: describes the microbial number parameter, either TVC, LAB (=lab) or fungi (=fy) counts  
**value:** TVC, LAB and fungi count  

# physical_chemical_temp.R

### data = physico_chemical.xlsx

**day:** measurment day. 0, 3, 6, 9 or 12. 0 marks the beginning of the feeding experiment, 12 the end of the feeding experiment     
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
