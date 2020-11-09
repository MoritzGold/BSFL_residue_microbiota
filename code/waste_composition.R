# Nutrient composition of rearing substrates (Table 1)

rm(list=ls()) 

library(extrafont)
library(openxlsx)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)
library(ggpubr)
library(readr)
library(knitr)
library(kableExtra)


#--------# 


# Read in data

  # all parameters except for non-fibre carbohydrates

  waste_messy <- read.xlsx(xlsxFile = "data/waste_composition.xlsx")
  
  # raw data for non-fibre carbohydrates (analyses performed with Megazyme enzyme kits)

  fructose <- read.xlsx(xlsxFile = "data/non_fibre_carbohydrates.xlsx",sheet = 1)
  glucose <- read.xlsx(xlsxFile = "data/non_fibre_carbohydrates.xlsx",sheet = 2)
  starch <- read.xlsx(xlsxFile = "data/non_fibre_carbohydrates.xlsx",sheet = 3)
  
  
# Determine non-fibre carbohydrate contents by conversion of absorbance values into % DM
  
  # fructose
  
    # determine absorbance difference of blanks
  
    fructose_blanks <-
    
    fructose %>% 
      
      filter(sample=="blank") %>% 
      
       mutate(blank_delta_A3_A2=A3-A2) %>% 
      
       select(blank_ref,blank_delta_A3_A2) 
  
    
    # determine absorbance differences for samples and substract absorbance of blanks
    
    fructose_perc_DM <- 
    
    fructose %>% 
      
      filter(!sample=="blank") %>% 
      
      mutate(delta_A3_A2=A3-A2) %>% 
      
      left_join(fructose_blanks,by="blank_ref") %>% 
      
      mutate(delta_A3_A2=delta_A3_A2-blank_delta_A3_A2) %>% 
      
      select(sample,Replicate,sample_weight_g,final_volume_mL,delta_A3_A2,moisture_content_perc) %>% 
      
      group_by(Replicate) %>% 
      
      # convert absorbance values into g/L
      # we use different treshhold for the absorbance difference, 0.65 instead of one
      
      mutate(fructose_g_L=case_when(delta_A3_A2>0.065 ~ 18.159*delta_A3_A2, TRUE~0)) %>% 
      
      # convert fructose in g/L into perc - see K-ACHDF Megazyme assay kit booklet for equation
      
      mutate(fructose_perc=case_when(fructose_g_L>0 ~ ((fructose_g_L/(sample_weight_g/(final_volume_mL/1000)))*100),
                                     TRUE~fructose_g_L)) %>% 
      
      # correct for moisture content in freeze-dried samples
      
      mutate(fructose_per_DM=fructose_perc*(1+(moisture_content_perc/100))) %>% 
      
      select(sample,Replicate,fructose_per_DM)
    
    
    # Calculate mean and difference between results
    
    fructose_perc_DM_diff <-
    
    fructose_perc_DM %>% 
      
      group_by(sample, Replicate) %>% 
      
      summarise(fructose_per_DM=mean(fructose_per_DM)) %>% 
      
      mutate(diff=lead(fructose_per_DM)) %>% 
      
      mutate(diff=diff-fructose_per_DM)
    
    
    fructose_perc_DM <- 
    
    fructose_perc_DM_diff %>% 
      
      mutate(value=fructose_per_DM) %>% 
      
      select(-diff,-fructose_per_DM) %>% 
      
      mutate(parameter="fructose")
    
    
    fructose_perc_DM %>% 
      
      summarise(mean=mean(value))
  
    
  # starch
    
    # determine absorbance difference of blanks
    
    starch_blanks <-
      
    starch %>% 
      
      filter(sample=="Control") %>% 
      
      group_by(sample,blank_ref) %>% 
      
      mutate(delta_A_blank=(delta_A_1+delta_A_2)/2) %>% 
      
      ungroup() %>% 
      
      select(blank_ref,delta_A_blank)
     
    
    # determine mean absorbance differences for samples
    
    starch_perc_DM <-
    
    starch %>% 
      
      left_join(starch_blanks,by=c("blank_ref")) %>% 
      
      mutate(delta_A=(delta_A_1+delta_A_2)/2) %>% 
      
      select(-blank_ref,-delta_A_1,-delta_A_2) %>% 
      
      # calculate F = absorbance of known glucose concentration
      
      mutate(F=100/delta_A_blank) %>% 
      
      filter(!sample=="Control") %>% 
      
      # calculate starch perc (final volume = 100 mL) - see K-TSTA Megazyme assay kit booklet for equation
      
      mutate(starch_perc=delta_A*(F/weight_mg)*100*0.9) %>% 
      
      # correct for moisture content in freeze-dried samples
      
      mutate(starch_per_DM=starch_perc*(1+(moisture_content_perc/100)))
  
    
   # Calculate mean and sd between results
    
    starch_perc_DM <-
    
    starch_perc_DM %>% 
      
      select(sample,Replicate,starch_per_DM) %>% 
      
      mutate(value=starch_per_DM) %>% 
      
      select(-starch_per_DM) %>% 
      
      mutate(parameter="starch")
    
    
    starch_perc_DM %>% 
      
      group_by(sample) %>% 
      
      summarise(mean=mean(value),
                sd=sd(value)) 
    
    
  # glucose
    
    # determine absorbance difference of blanks
    
    glucose_blanks <-
      
      glucose %>% 
      
      filter(sample=="Control") %>% 
      
      group_by(sample,blank_ref) %>% 
      
      mutate(delta_A_blank=(delta_A_1+delta_A_2)/2) %>% 
      
      ungroup() %>% 
      
      select(blank_ref,delta_A_blank)
      
    
  # determine mean absorbance differences for samples
    
    glucose_perc_DM <-
      
      glucose %>% 
      
      left_join(glucose_blanks,by=c("blank_ref")) %>% 
      
      mutate(delta_A=(delta_A_1+delta_A_2)/2) %>% 
      
      select(-blank_ref,-delta_A_1,-delta_A_2) %>% 
        
      filter(!sample=="Control") %>%   
      
      # calculate glucose ug / 0.1 mL - see K-GLUC Megazyme assay kit booklet for equation 
        
      mutate(glucose_ug_0.1mL=(delta_A/delta_A_blank)*100) %>% 
        
      # convert glucose ug / 0.1 mL into perc
        
      mutate(glucose_perc=((glucose_ug_0.1mL*volume_mL*10)/weight_mg)/10)  %>% 
      
      # correct for moisture content in freeze-dried samples
      
      mutate(glucose_perc_DM=glucose_perc*(1+(moisture_content_perc/100))) %>% 
      
      select(sample,Replicate,glucose_perc_DM)

    
   # Calculate mean and sd between results
    
    glucose_perc_DM <-
      
      glucose_perc_DM %>% 
      
      mutate(parameter="glucose") %>% 
      
      mutate(value=glucose_perc_DM) %>% 
      
      select(-glucose_perc_DM)
    
    
    glucose_perc_DM %>% 
      
      group_by(sample) %>% 
      
      summarise(mean=mean(value),
                sd=sd(value))  
        

    
  # non-fibre carbohydrate summary
    
    NFC <- 
    
    bind_rows(fructose_perc_DM,
              glucose_perc_DM,
              starch_perc_DM) %>% 
      
      group_by(sample,Replicate) %>% 
      
      spread(parameter,value) %>% 
      
      na.omit() %>% 
      
      mutate(NFC=fructose+glucose+starch)
    
    
    NFC %>% 
      
      select(sample,NFC) %>% 
      
      ungroup() %>% 
      
      mutate(diff=lead(NFC)) %>% 
      
      mutate(diff=diff-NFC)
    
    
    NFC_mean <-
    
    NFC %>% 
      
      group_by(sample) %>% 
      
      summarise(mean=mean(NFC)) %>% 
      
      mutate(Parameter="NFC") %>% 
      
      rename(Diet=sample) 
 
           
#Convert measurements into nutrient compostion metrics

waste_composition <- 

waste_messy %>% 
  
  select(Diet,Replicate,Parameter,Value) %>% 
  
  spread(key = "Parameter",value = "Value") %>% 
  
  #Convert nitrogen into protein - we use factor of 5.4 (see Gold et al. (2020), Waste Management, 102, 319-329)
  
  mutate(Protein=Nitrogen*5.4) %>% 
  
  select(-Nitrogen) %>% 
  
  #Estimate hemicelluloses by difference between NDF and ADF
  
  mutate(Hemicellulose=NDF-ADF) %>% 
  
  select(-NDF) %>% 
  
  mutate(Cellulose_lignin=ADF) %>% 
  
  select(-ADF,-Ash,-CF) %>% 
  
  #Calculate fibre by sum of hemicellulose and cellulose and lignin
  
  mutate(Fibre=Cellulose_lignin+Hemicellulose) %>% 

  gather(3:8,key = Parameter,value = Value,na.rm = TRUE)
  
  


#calculate means and show replicate measurments

waste_composition_mean <- 

waste_composition %>% 
  
  group_by(Diet,Parameter) %>% 
  
  summarise(n= n(), 
            mean=mean(Value))



#Calculate standard deviation (n>2) and difference between measurements (n=2)


  #standard deviation

  waste_composition %>% 
  
    left_join(waste_composition_mean,by=c("Diet","Parameter")) %>% 
    
    select(-mean) %>% 
    
    filter(n>2) %>%
    
    group_by(Diet,Parameter) %>% 
    
    summarise(sd=sd(Value))
    
    
  #difference between measurements
  
  waste_composition %>% 
    
    left_join(waste_composition_mean,by=c("Diet","Parameter")) %>% 
    
    select(-mean) %>% 
    
    filter(n==2) %>%
    
    mutate(lead=lead(Value)) %>% 
    
    mutate(diff=lead-Value) %>% 
    
    filter(Replicate==1) %>% 
    
    select(Diet,Parameter,diff)
    

  
#Estimate caloric content
  
  waste_composition_mean %>% 
    
    select(-n) %>% 
    
    bind_rows(NFC_mean) %>% 
    
    filter(Parameter=="Lipids"|Parameter=="NFC"|Parameter=="Protein") %>% 
    
    group_by(Diet) %>% 
    
    spread(3,key = Parameter,value = mean) %>% 
    
    mutate(caloric_content=(9.4*Lipids)+(5.4*Protein)+(4.1*NFC))
