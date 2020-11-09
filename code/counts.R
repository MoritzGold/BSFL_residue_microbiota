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

counts <- read.xlsx(xlsxFile = "data/counts.xlsx")


#unit CFU/g residue


#residue


# calculate mean and sd (n>2)

  counts_residue <-
    
    counts %>% 
    
    filter(type=="residue") %>% 
    
    group_by(day, type, biowaste,parameter) %>% 
    
    mutate(value=log(value,base = 10))


  counts_residue_mean <-
    
  counts %>% 
    
    filter(type=="residue") %>% 
    
    group_by(day, type, biowaste,parameter) %>% 
    
    mutate(value=log(value,base = 10)) %>% 
    
    summarise(n=n(),
              mean=mean(value),
              sd=sd(value)) %>% 
    
  mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_))
  

# saved data to MS Excel

  counts_residue_mean_write <-    
    
    counts_residue_mean %>% 
    
    gather(5:7,key = stat,value=value) %>% 
    
    spread(key = parameter,value = value)
  
  
  write.xlsx(counts_residue_mean_write, file="output/residue_counts.xlsx")


# calculate diff (n=2)

  counts_residue_diff <-
  
  counts_residue_mean %>% 
    
    filter(n==2) %>% 
    
    select(day,type,biowaste,parameter) %>% 
    
    left_join(counts_residue,by=c("day","type","biowaste","parameter")) %>% 
    
    mutate(diff=lag(value)) %>% 
    
    mutate(diff=value-diff)

  
  write.xlsx(counts_residue_diff, file="output/residue_counts_diff.xlsx")
  
#plot

counts_residue_mean %>% 

  ggplot(aes(day,mean)) +
  
  geom_point(aes(shape=biowaste),size=3) +
  
  geom_line(aes(group=biowaste),size=.2,linetype="dashed") + 
  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)+
  
  geom_jitter(data=counts_residue,
              aes(day,value),
              size=2,
              position = position_jitter(0.2),
              alpha=.4) +
  
  theme_bw(base_size = 18)+
  
  facet_wrap(~parameter) +
  
  ylim(0,12)
