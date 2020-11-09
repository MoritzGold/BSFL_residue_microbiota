#Physical-chemical data and temperature of residue

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


#--------# Read in data

physchem_messy <- read.xlsx(xlsxFile = "data/physico_chemical.xlsx")

temp_messy <- read.xlsx(xlsxFile = "data/temp.xlsx")



#--------# Calculate interesting parameters

psychem <- 

physchem_messy %>%

  group_by(day,biowaste,type,replicate,parameter) %>% 
  
  summarise(value=mean(value,na.rm = TRUE)) %>% 
  
  spread(key = parameter,value = value) %>% 
  
  mutate(C_N=C/N) %>% 
  
  mutate(VS_ash=VS/ash) %>% 
  
  gather(5:13,key = parameter,value = value,na.rm = TRUE)



#--------# write to table
 
  # mean and sd (n>2)

  psychem_table <-
  
    psychem %>%
    
    group_by(day,biowaste,type,parameter) %>% 
    
    summarise(n=n(),
              mean=mean(value,na.rm=TRUE),
              sd=sd(value,na.rm=TRUE)) %>% 
      
    mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_)) %>% 
     
    filter(parameter=="aw"|parameter=="C_N"|parameter=="MC"|parameter=="N"|parameter=="pH"|parameter=="C"|parameter=="VS")
  
  write.xlsx(psychem_table, file="output/psychem.xlsx")


  # diff (n=2)
  
  psychem_table_diff <- 
  
  psychem_table %>% 
    
    filter(n==2) %>% 
    
    select(day,biowaste,type,parameter) %>% 
    
    left_join(psychem,by=c("day", "biowaste","type","parameter")) %>% 
    
    mutate(diff=lag(value)) %>% 
    
    mutate(diff=value-diff)

    
  write.xlsx(psychem_table_diff, file="output/psychem_diff.xlsx")

  
#--------# temperature


temp_mean <- 

temp_messy %>% 
  
  group_by(day,biowaste,type) %>% 
  
  summarise(n=n(),
            mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE)) %>% 
  
  mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_)) %>% 
  
  select(-n)


 write.xlsx(temp_mean, file="output/temp_residue.xlsx")

 
 temp_mean %>% 
   
   ungroup() %>% 
   
   summarise(min=min(mean),
             max=max(mean))

  #--------# Correlation between parameters 


# prepare data for use in phyloseq.R

psychem_C_H <- 

psychem %>% 
  
  filter(biowaste=="C"|biowaste=="H") %>% 
  
  filter(!type=="no_larvae") %>% 
  
  spread(key = parameter,value = value) %>% 
  
  ungroup() %>% 
  
  select(-type) %>% 
  
  rename(Day = day) %>% 
  
  rename(Diet = biowaste) %>% 
  
  rename(Replicate = replicate) %>% 
  
  select(-ash,-VS_ash)
  
 

# correlations between rearing performance and physio-chemical residue compostion 

  # see Pearson correlation plot in phyloseq.R, used to identify co-linearity
  # between parameters for dbRDA. 


  # confirm data is similiarly distributed to a normal distribution

  jpeg("output/larval_weight_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$DM_Larvae_mg, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$DM_Larvae_mg, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n Larval weight mg DM", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()
  
  jpeg("output/WR_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$Waste_reduction_DM, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$Waste_reduction_DM, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n Waste reduction % DM", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()
  
  jpeg("output/BCR_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$bioconversion_rate, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$bioconversion_rate, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n Bioconversion rate % DM", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()
  
  jpeg("output/protein_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$P_gDM_container, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$P_gDM_container, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n Protein content g DM/replicate", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()
  
  jpeg("output/MC_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$MC, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$MC, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n Moisture content % DM", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()
  
  jpeg("output/pH_corr.jpeg", width = 350, height = 350,quality = 1200)
  qqnorm(performance_psychem$pH, pch = 1, frame = FALSE,main = NA)
  qqline(performance_psychem$pH, col = "steelblue", lwd = 2)
  mtext("Q-Q plot\n pH", cex = 1.2,side = 3, line = -2.5, outer = TRUE)
  dev.off()


  # residue moisture content
  
  cor.test(performance_psychem$DM_Larvae_mg, 
           performance_psychem$MC,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$Waste_reduction_DM, 
           performance_psychem$MC,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$bioconversion_rate, 
           performance_psychem$MC,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$P_gDM_container, 
           performance_psychem$MC,  method = "pearson",exact = FALSE)


  # residue pH
  
  cor.test(performance_psychem$DM_Larvae_mg, 
           performance_psychem$pH,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$Waste_reduction_DM, 
           performance_psychem$pH,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$bioconversion_rate, 
           performance_psychem$pH,  method = "pearson",exact = FALSE)
  
  cor.test(performance_psychem$P_gDM_container, 
           performance_psychem$pH,  method = "pearson",exact = FALSE)
  
  

  