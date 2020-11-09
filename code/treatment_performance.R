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
library(cowplot)
library(gridExtra)
library(vegan)


#--------# 


#Read in data ----

performance_messy <- read.xlsx(xlsxFile = "data/treatment_performance.xlsx")



#Larval weight ----

larval_weight_mean <- 

performance_messy %>% 
  
  select(Day,Type,Diet,Replicate,DM_Larvae_mg) %>% 
  
  group_by(Diet,Day,Type) %>% 
  
  summarise(n=n(),
            mean=mean(DM_Larvae_mg),
            sd=sd(DM_Larvae_mg)) %>% 
  
  mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_))


larval_weight <-

  performance_messy %>% 
  
  filter(!(Type=="Control C"|Type=="Control S")) %>% 
  
  select(Day,Type,Diet,Replicate,DM_Larvae_mg)



#Waste reduction ----

wr_mean <- 
  
  performance_messy %>% 
  
  select(Day,Diet,Type,Replicate,Waste_reduction_DM) %>% 
  
  group_by(Diet,Day,Type) %>% 
  
  summarise(n=n(),
            mean=mean(Waste_reduction_DM),
            sd=sd(Waste_reduction_DM)) %>% 
  
  mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_))


wr <-
  
  performance_messy %>% 
  
  filter(!(Type=="Control C"|Type=="Control S")) %>% 
  
  select(Day,Type,Diet,Replicate,Waste_reduction_DM)


##Bioconversion rate ----
  
  
  BCR <-
  
    performance_messy %>% 
    
    filter(!(Type=="Control C"|Type=="Control S")) %>% 
    
    select(Day,Diet,Replicate,Number_larvae_in,DM_Larvae_mg,Feed_in_gDM) %>% 
    
    mutate(bioconversion_rate = (((((Number_larvae_in*DM_Larvae_mg)-(Number_larvae_in*0.5))/1000)/Feed_in_gDM)*100))
  
  
  BCR_mean <- 
  
  performance_messy %>% 
    
  select(Day,Diet,Replicate,Type,Number_larvae_in,DM_Larvae_mg,Feed_in_gDM) %>% 
    
  mutate(bioconversion_rate = (((((Number_larvae_in*DM_Larvae_mg)-(Number_larvae_in*0.5))/1000)/Feed_in_gDM)*100)) %>% 
    
  group_by(Diet,Day,Type) %>% 
    
  summarise(n=n(),
              mean_bcr=mean(bioconversion_rate),
              sd=sd(bioconversion_rate)) %>% 
    
    mutate(sd=case_when(n > 2 ~ sd, TRUE~ NA_real_))
 
  

#Larval composition ----
  
  
  # mean larval composition
  
  larval_composition_mean <- 
    
    performance_messy %>% 
    
    filter(Type=="Sample") %>% 
    
    select(Day,Diet,Replicate,Ash_Larvae_percDM,DM_Larvae_mg,N_Larvae_percDM) %>% 
    
    mutate(P_Larvae_percDM=N_Larvae_percDM*4.67) %>% 
    
    select(-N_Larvae_percDM) %>% 
    
    mutate(P_gDM_container=((200*((P_Larvae_percDM/100)*DM_Larvae_mg))/1000)) %>%
    
    select(-DM_Larvae_mg) %>% 
    
    gather(4:6,key = parameter,value = value) %>% 
    
    group_by(Diet,Day,parameter) %>% 
    
    summarise(n=n(),
              mean_comp=mean(value),
              sd_comp=sd(value)) %>% 
    
    mutate(sd_comp=case_when(n > 2 ~ sd_comp, TRUE~ NA_real_))
  
  
  # larval compostion per replicate
  
  larval_composition <-
    
    performance_messy %>% 
    
    filter(Type=="Sample") %>% 
    
    select(Day,Diet,Replicate,Ash_Larvae_percDM,DM_Larvae_mg,N_Larvae_percDM) %>% 
    
    mutate(P_Larvae_percDM=N_Larvae_percDM*4.67) %>% 
    
    select(-N_Larvae_percDM) %>% 
    
    mutate(P_gDM_container=((200*((P_Larvae_percDM/100)*DM_Larvae_mg))/1000)) %>%
    
    select(-DM_Larvae_mg) %>% 
    
    gather(4:6,key = parameter,value = value)
    
 
  # mean larval protein accumulation
  
  larval_protein_accumulation_mean <- 
  
  larval_composition_mean %>% 
    
    filter(parameter=="P_gDM_container")
    
  
  # larval protein accumulation per replicate
  
  larval_protein_accumulation <- 
  
  larval_composition %>% 
    
    filter(parameter=="P_gDM_container") %>% 
    
    spread(key = parameter,value = value)

  
  
#Data wrangling for stacked bar graphs ----
 
  
  #Treatments without BSFL ----
  
  performance_sum_high_n_control <-
    
    performance_messy %>% 
    
    select(Day,Diet,Type,Replicate,Waste_reduction_DM) %>% 
    
    filter((Type=="Control C"|Type=="Control S")) %>% 
    
    group_by(Type) %>% 
    
    summarise(n=n(),
              sd=sd(Waste_reduction_DM),
              mean=round(mean(Waste_reduction_DM),1)) %>%    
    
    mutate(n=case_when(n>2 ~ "",
                       n==2 ~ "†",
                       TRUE~"‡")) %>% 
                       
    mutate(value_new=mean) %>% 
    
    unite(label,value_new,n,sep=", ") %>% 
    
    mutate(value=mean) %>% 
    
    mutate(Day="day: 0-12") %>% 
    
    mutate(parameter="WR_percDM") %>% 
    
    mutate(Diet=case_when(Type=="Control C"~"FW-0",
                          Type=="Control S"~"S-FW-0",
                          TRUE~"error")) %>% 
    
    select(-Type)
  
  

  #Treatment with BSFL ----
   

  #summarise mean process performance values

  performance_sum <-
  
  BCR_mean %>% 
    
    filter(!(Type=="Control C"|Type=="Control S")) %>% 
    
    select(-n,-sd)%>% 
    
    left_join(larval_weight_mean,by=c("Diet","Day","Type")) %>% 
    
    select(-n,-sd) %>% 
    
    left_join(wr_mean,by=c("Diet","Day","Type")) %>% 
    
    select(-n,-sd,-Type) %>% 
    
    left_join(larval_protein_accumulation_mean,by=c("Diet","Day")) %>% 
    
    select(-n,-sd_comp,-parameter)
  

  #change column names  

  names(performance_sum) = c("Diet","Day","BCR_percDM","larval_mgDM","WR_percDM","P_gDM_container")  
    
  
  #gather and round values
  
  performance_mean <- 

  performance_sum %>% 
    
    gather(3:6,key = parameter,value = mean) %>% 
    
    na.omit(mean)
    
  
  #calculate differences in performance between measurements days
  
  
    #for performance metrics without a value at the beginning of the rearinge experiment (larval DM and protein content)
  
    performance_sum_high_P_larvae_DM <- 
         
      performance_sum %>% 
      
      gather(3:6,key = parameter,value = value) %>% 
  
      filter(parameter=="P_gDM_container"|parameter=="larval_mgDM") %>% 
          
      group_by(Diet, parameter) %>% 
      
      mutate(lag=lag(value,order_by = Day,n = 1)) %>% 
      
      mutate(value_diff=value-lag) %>% 
      
      mutate(value_diff=case_when(Day=="0"~value,TRUE~value_diff)) %>% 

      select(-value, -lag) %>% 
      
      rename(value = value_diff)
      
 
    #for performance metrics without a value at the beginning of the rearinge experiment (WR and BCR)
       
    performance_sum_high_WR_BCR <- 
      
    performance_sum %>% 
      
      gather(3:6,key = parameter,value = value) %>% 
      
      group_by(Diet, parameter) %>% 
  
      filter(!parameter=="P_gDM_container"&!parameter=="larval_mgDM") %>% 
      
      mutate(lag=lag(value,order_by = Day,n = 1)) %>% 
      
      mutate(value_diff=value-lag) %>% 
      
      mutate(value_diff=case_when(Day=="3"~value,TRUE~value_diff)) %>% 
      
      filter(Day>0) %>% 
  
      select(-value, -lag) %>% 
      
      rename(value = value_diff)
    
    
    #merge tables
      
    performance_sum_high <- bind_rows(performance_sum_high_WR_BCR, performance_sum_high_P_larvae_DM) 
  

  #summarise number of replicates and sd

  BCR_n_sd <- 
   
  BCR_mean %>% 
    
    filter(!(Type=="Control C"|Type=="Control S")) %>% 
    
    select(-Type) %>% 
    
    filter(Day>0) %>% 
    
    select(-mean_bcr) %>% 
    
    mutate(parameter="BCR_percDM")
  
  
  larva_n_sd <- 
    
    larval_weight_mean %>% 
    
    filter(!(Type=="Control C"|Type=="Control S")) %>% 
    
    select(-Type) %>% 
    
    select(-mean) %>% 
    
    mutate(parameter="larval_mgDM")
  
  
  wr_n_sd <-
  
  wr_mean %>% 
    
    filter(Day>0) %>% 
    
    filter(!(Type=="Control C"|Type=="Control S")) %>% 
    
    select(-Type) %>% 
    
    select(-mean) %>% 
    
    mutate(parameter="WR_percDM")
  
  
  P_n_sd <-
    
    larval_protein_accumulation_mean %>% 
    
    select(-mean_comp) %>% 
    
    rename(sd = sd_comp)

  
  sum_s_n <- bind_rows(BCR_n_sd,larva_n_sd, wr_n_sd, P_n_sd) %>% 
    
    #convert n to icons to illustrate number of replicates
    
    mutate(n_new=n) %>% 
      
    mutate(n=case_when(n>2 ~ "",
                       n==1 ~ "†",
                       n==2 ~ "‡",
                       TRUE~"error")) 

  #add sd and n to the matrix including the change in performance metrics between measurement days
  
  performance_sum_high <-
  
  performance_sum_high %>% 
    
    left_join(sum_s_n, by=c("Diet","Day","parameter")) %>% 
  
    mutate(value_new=round(value,1)) %>% 
    
    unite(label,value_new,n,sep="")
  
  
  #changing the diet and day variables and adding mean value of performance metrics

  #day_range     <- c("day: 0-3","day: 4-6","day: 7-9","day: 10-12")
  day_range     <- c("day < 0","day: 0-3","day: 4-6","day: 7-9","day: 10-12","day: 0-12")


  performance_sum_high_n <-
  
  performance_sum_high %>% 
    
    left_join(performance_mean, by=c("Diet","Day","parameter")) %>% 
    
    mutate(Day=case_when(Day=="0"~"day < 0",
                         Day=="3"~"day: 0-3",
                         Day=="6"~"day: 4-6",
                         Day=="9"~"day: 7-9",
                         Day=="12"~"day: 10-12",TRUE~"error")) %>% 

    #bind with results of treatments without BSFL
  
    bind_rows(performance_sum_high_n_control) %>% 
    
    mutate(Day=factor(Day,levels=day_range)) %>% 
    
    ungroup(Diet) %>% 
    
    mutate(Diet=case_when(Diet=="C"~"FW",
                          Diet=="S"~"S-FW",
                          Diet=="H"~"HW",TRUE~Diet))
 
  
  
#Plot stacked bar graphs ---- 
    
  
  stacked_bar_larva <-
    
    performance_sum_high_n %>% 
    
    mutate(Diet=case_when(Diet=="FW"~"Canteen\n waste",
                          Diet=="S-FW"~"Sterile\ncanteen\nwaste",
                          Diet=="HW"~"Household\nwaste",TRUE~Diet)) %>% 
    
    filter(parameter=="larval_mgDM") %>% 
    
    ggplot(aes(x = Diet, y = value, fill = Day,label=label)) +
    
    geom_bar(position = position_stack(reverse = TRUE), stat="identity") +
    
    geom_bar(stat = "identity",position = position_stack(reverse = TRUE),color="white") +
    
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2,alpha=.4) +
    
    geom_text(size = 3, position = position_stack(vjust = 0.5,reverse = TRUE)) +
    
    theme_bw(base_size = 14) +
    
    xlab("") +
    
    ylab("Larval weight\n[mg DM / larva]") +
    
    scale_fill_manual(name="Rearing\nduration",values = c("#661100","#888888","#117733","#CC6677","#88CCEE")) +
    
    labs(tag = "A")



stacked_bar_protein <-

performance_sum_high_n %>% 
  
  mutate(Diet=case_when(Diet=="FW"~"Canteen\n waste",
                        Diet=="S-FW"~"Sterile\ncanteen\nwaste",
                        Diet=="HW"~"Household\nwaste",TRUE~Diet)) %>% 
  
  filter(parameter=="P_gDM_container") %>% 

  ggplot(aes(x = Diet, y = value, fill = Day,label=label)) +
  
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE),color="white") +
  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2,alpha=.4) +
  
  geom_text(size = 3, position = position_stack(vjust = 0.5,reverse = TRUE)) +
  
  theme_bw(base_size = 14) +
  
  xlab("") +
  
  ylab("Larval protein\n[g DM / replicate]") +
  
  scale_fill_manual(name="Rearing\nduration",values = c("#661100","#888888","#117733","#CC6677","#88CCEE")) +
  
  labs(tag = "C") +
  
  theme(legend.position="none")
 

stacked_bar_BCR <-

performance_sum_high_n %>% 
  
  filter(parameter=="BCR_percDM") %>% 
  
  mutate(Diet=case_when(Diet=="FW"~"Canteen\n waste",
                        Diet=="S-FW"~"Sterile\ncanteen\nwaste",
                        Diet=="HW"~"Household\nwaste",TRUE~Diet)) %>% 
  
  ggplot(aes(x = Diet, y = value, fill = Day,label=label)) +
  
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE),color="white") +
  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2,alpha=.4) +                
  
  geom_text(size = 3, position = position_stack(vjust = 0.5,reverse = TRUE)) +
  
  theme_bw(base_size = 14) +
  
  ylab("Bioconversion rate\n[% DM]") +
  
  xlab("") +
  
  scale_fill_manual(name="Rearing\nduration",values = c("#888888","#117733","#CC6677","#88CCEE")) +
  
  labs(tag = "B") +
  
  theme(legend.position="none")


stacked_bar_WR_1 <-

performance_sum_high_n %>% 
  
  filter(Diet=="FW"|Diet=="HW") %>% 
  
  mutate(Diet=case_when(Diet=="FW"~"Canteen\n waste\n\n",
                        Diet=="HW"~"Household\nwaste\n\n",TRUE~Diet)) %>% 
  
  filter(parameter=="WR_percDM") %>% 
  
  ggplot(aes(x = Diet, y = value, fill = Day,label=label)) +
  
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE),color="white") +
  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2,alpha=.4) +
  
  geom_text(size = 3, position = position_stack(vjust = 0.5,reverse = TRUE)) +
  
  theme_bw(base_size = 14) +
  
  ylab("Waste reduction\n[% DM]") +
  
  xlab("") +
  
  scale_fill_manual(name="Rearing\nduration",values = c("#888888","#117733","#CC6677","#88CCEE")) +
  
  expand_limits(y = c(-5, 70)) +
  
  labs(tag = "A") +
  
  theme(legend.position="none")



stacked_bar_WR_2 <-
  
  performance_sum_high_n %>% 
  
  filter(parameter=="WR_percDM") %>% 
  
  filter(Diet=="FW"|Diet=="S-FW"|Diet=="FW-0"|Diet=="S-FW-0") %>% 
  
  mutate(Diet=case_when(Diet=="FW"~"Canteen\n waste",
                        Diet=="HW"~"Household\nwaste",
                        Diet=="S-FW"~"Sterile\ncanteen\nwaste",
                        Diet=="FW-0"~"Canteen\nwaste\n(without BSFL)",
                        Diet=="S-FW-0"~"Sterile\ncanteen\nwaste\n(without BSFL)",
                        TRUE~Diet)) %>% 
  
  ggplot(aes(x = Diet, y = value, fill = Day,label=label)) +
  
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE),color="white") +
  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.2,alpha=.4) +
  
  geom_text(size = 3, position = position_stack(vjust = 0.5,reverse = TRUE)) +
  
  theme_bw(base_size = 14) +
  
  ylab("Waste reduction\n[% DM]") +
  
  xlab("") +
  
  scale_fill_manual(name="Rearing\nduration",values = c("#888888","#117733","#CC6677","#88CCEE","#DDCC77")) +
  
  labs(tag = "B") +
  
  expand_limits(y = c(-5, 70))

#combine the individual plots


#get legend

legend_1 <- get_legend(stacked_bar_larva)
legend_2 <- get_legend(stacked_bar_WR_2)
                       
#remove legend

stacked_bar_larva <- stacked_bar_larva + theme(legend.position="none")
stacked_bar_WR_2 <- stacked_bar_WR_2 + theme(legend.position="none")

#combine plots and produce grid

bar_grid_1 <-  grid.arrange(stacked_bar_larva,stacked_bar_BCR,stacked_bar_protein,legend_1,ncol=4,nrow=1,widths = c(3,3,3,1))

bar_grid_2 <-  grid.arrange(stacked_bar_WR_1,stacked_bar_WR_2,legend_2,ncol=3, nrow=1,widths = c(3.5, 5.5,1.5))
                            

save_plot("output/performance_bar_1.jpeg", bar_grid_1,base_height = 8, base_width=12)
save_plot("output/performance_bar_2.jpeg", bar_grid_2,base_height = 8, base_width=9)





#--------# Correlation between parameters ----


# prepare data for use in phyloseq.R


performance_C_H <- 
  
  performance_messy %>% 
  
  filter(Diet=="C"|Diet=="H") %>% 
  
  filter(!Type=="Control C") %>% 
  
  select(2,4,5,8:12,15) %>% 
  
  left_join(BCR,by=c("Day","Diet","Replicate","DM_Larvae_mg")) %>% 
  
  left_join(larval_protein_accumulation, by=c("Day","Diet","Replicate")) %>% 

  select(-N_Larvae_percDM,-C_Larvae_percDM,-Ash_Larvae_percDM,-VS_Larvae_percDM,-Number_larvae_in,
       -Feed_in_gDM)



