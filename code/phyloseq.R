# Libraries and functions -------------------------------------------------

## clean/reset environment 

  rm(list=ls()) 

## R and Bioconductor libraries 

  require("ggplot2")
  require("vegan")
  #require("knitr")
  require("phyloseq")
  require("microbiome")
  #require("ampvis2")
  #require("betapart")
  #require("dplyr")
  #require("magrittr")
  #require("doBy")
  #require("plotly")
  #require("cooccur")
  #require("scales")
  #require("picante")
  library(tidyverse)
  #library(metacoder)
  library(ampvis2)
  library(openxlsx)
  library(UpSetR)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggord)
  library(corrplot)
  library(usdm)
  library(conflicted)
  library(gridExtra)  

  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("intersect", "dplyr")
  conflict_prefer("lag", "dplyr")

  
## Functions

  source("code/Rfunctions.R")

# Data import in R --------------------------------------------------------

  zotufile    <- "code/p529_run190225_16S_ZOTU_Count_Sintax.txt"
  mapfile     <- "code/p529_run190225_MapFile.txt"
  treefile    <- "code/p529_run190225_16S_ZOTU_MSA.tre"
  
  d <- import_qiime(otufilename = zotufile, mapfilename = mapfile, treefilename = treefile)

# Remove rare taxa < 10 reads ---------------------------------------------------------
  
  d <- prune_taxa(taxa_sums(d) >= 10, d)    
  
  
# Convert reads in relative abundance per sample ----------------------------------------------
  
  d_norm <- transform_sample_counts(d, function(x) x / sum(x))

  
# Mock culture -----
  
  
  ### Expected Diversity of the Mock
  
  # 1  Bacteria	Actinobacteria	Actinobacteria	    Micrococcales	    Micrococcaceae	    Arthrobacter	Pseudarthrobacter chlorophenolicus
  # 2  Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales  	Burkholderiaceae	  Burkholderia	Burkholderia xenovorans
  # 3  Bacteria	Firmicutes	    Bacilli            	Bacillales	      Bacillaceae	        Bacillus	    Bacillus subtilis
  # 4  Bacteria	Proteobacteria	Zammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	  Escherichia coli
  # 5  Bacteria	Actinobacteria	Actinobacteria    	Micrococcales	    Micrococcaceae	    Micrococcus	  Micrococcus luteus
  # 6  Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales  	Pseudomonadaceae   	Pseudomonas	  Pseudomonas protegens
  # 7  Bacteria	Firmicutes	    Bacilli	            Bacillales	      Paenibacillaceae	  Paenibacillus	Paenibacillus sabinae
  # 8  Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	  Pseudomonadaceae  	Pseudomonas	  Pseudomonas stutzeri
  # 9  Bacteria	Actinobacteria	Actinobacteria	    Streptomycetales	Streptomycetaceae	  Streptomyces	Streptomyces violaceoruber
  # 10 Bacteria	Proteobacteria	Alphaproteobacteria	Rhizobiales	      Xanthobacteraceae	  Xanthobacter	Xanthobacter autotrophicus
  
  d_mock <- subset_samples(d, Type == "Mock")
  
  sort(taxa_sums(d_mock), TRUE) 
  
  #ZOTU32  ZOTU30  ZOTU57  ZOTU58  ZOTU60  ZOTU62  ZOTU87 ZOTU117 ZOTU116 ZOTU154  
  #10821   10114    9349    7861    7571    7484    2973    1238    1231     646
  
  tax_table(d_mock)[c("ZOTU32","ZOTU30","ZOTU57","ZOTU58","ZOTU60","ZOTU62","ZOTU87","ZOTU117","ZOTU116","ZOTU154")]
  
  #3 #ZOTU32  "Firmicutes"     "Bacilli"             "Bacillales"        "Bacillaceae"        "Bacillus"                         
  #4 #ZOTU30  "Proteobacteria" "Gammaproteobacteria" "Enterobacteriales" "Enterobacteriaceae" "Escherichia-Shigella"               
  #6/8 #ZOTU57  "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales"  "Pseudomonadaceae"   "Pseudomonas"                       
  #5(1) #ZOTU58  "Actinobacteria" "Actinobacteria"      "Micrococcales"     "Micrococcaceae"     "Micrococcus"                        
  #9 #ZOTU60  "Actinobacteria" "Actinobacteria"      "Streptomycetales"  "Streptomycetaceae"  "Streptomyces"                      
  #2 #ZOTU62  "Proteobacteria" "Betaproteobacteria"  "Burkholderiales"   "Burkholderiaceae"   "Burkholderia-Paraburkholderia"     
  #7 #ZOTU87  "Firmicutes"     "Bacilli"             "Bacillales"        "Paenibacillaceae"   "Paenibacillus"                     
  #7 #ZOTU117  Firmicutes"     "Bacilli"             "Bacillales"        "Paenibacillaceae"   "Paenibacillus"                      
  #3 #ZOTU116 "Firmicutes"     "Bacilli"             "Bacillales"        "Bacillaceae"        "Bacillus"                           
  #7 #ZOTU154  "Firmicutes"     "Bacilli"             "Bacillales"        "Paenibacillaceae"   "Paenibacillus"                      
  
  # Missing reference:
  # 10 Bacteria	Proteobacteria	Alphaproteobacteria	Rhizobiales	      Xanthobacteraceae	  Xanthobacter	Xanthobacter autotrophicus
  
  
  d_mock_norm <- subset_samples(d_norm, Type == "Mock")
  
  d_mock_norm <- filter_taxa(d_mock_norm, function(x) mean(x) > 0.01, TRUE)
  
  otu_table(d_mock_norm)
  
  # ZOTU2 is abundant but should not be present in the Mock culture - look out for this taxa during result interpretation

  
# Look at controls -------------------------------------------------

  
  # This study included both controls in DNA extraction (Type=control_DNA) and during library preparation (control_lib_prep).
  
  d_controls <- subset_samples(d, Type=="control_DNA" | Type =="control_lib_prep")
  
  summary(sample_sums(d_controls)) 
  
  sample_sums(d_controls)
  
  # All control_DNA have a high number of reads. 
  # 2/3 of control_lib_prep have a high number of reads. 
  # This is problematic and needs to be further evaluated before data analysis.
  
  
#---
  
  
  # Look seperately at library prep controls processed next to larvae and residue samples
  
    #larvae
  
    d_controls_larvae <- subset_samples(d, Type_detail =="control_lib_prep_larvae")
    
    summary(sample_sums(d_controls_larvae)) 
    
    sample_sums(d_controls_larvae)

    #residue
    
    d_controls_residue <- subset_samples(d, Type_detail =="control_lib_prep_residue")
    
    summary(sample_sums(d_controls_residue)) 
    
    sample_sums(d_controls_residue)
  
  # All residue library prep controls have low number of reads
  # All larvae library prep controls have high number of reads

    

#---

         
     
  # Where does the contamination come from? Look at beta diversity (UniFrac distance, NMDS).    
    
    d_overview <-  subset_samples(d)
    
    day <- c("0","3","6","9","12")
    
    d_overview_new <- phyloseq::filter_taxa(d_overview,function(x) sum(x) > 0.01,TRUE)
    
    sample_data(d_overview)$Day <- factor(x=sample_data(d_overview)$Day, levels=day)
    
    get_variable(d_overview, "Day")
    
    d_overview_BC <- ordinate(d_overview, "NMDS", "unifrac", weighted=TRUE)
    
    plot_ordination(d_overview_new, d_overview_BC,type="samples",color="Type",
                    title="") + 
      
      geom_point(size=4,alpha=.5) +
      
      geom_text(mapping = aes(label = Day), size = 3, vjust = 1.5) +
      
      theme_bw(base_size = 15) +
      
      scale_color_discrete(name  ="Sample type",
                           breaks=c("control_DNA", "control_lib_prep","larvae","Mock","poulty_feed_residue","residue"),
                           labels=c("Control DNA extraction", "Control library prep","larvae","Mock","Poultry feed residue","Residue"))
    
  # Reads in controls seem to be associated with larval samples. 

    
#---
   
    
# Remove of samples ------------------- 
  
    
  # Larval samples - possible cross-contamination
  
  d <- subset_samples(d, Type == "residue" | Type=="control_DNA" | Type_detail =="control_lib_prep_residue")  


  # chicken feed substrate - not relevant
  
  d <- subset_samples(d, !Diet == "CF")

  
  # phylum cynaobacterium (=plant chloroplasts)
  
  d = subset_taxa(d, !Phylum=="Cyanobacteria")

  
  # family Mitochondria (=plant mitochondria)
  
  d = subset_taxa(d, !Family=="Mitochondria")

  
  # Residue, canteen waste, day 12, replicate 3 - deviates very much from the other three samples

  d <- subset_samples(d, !(Type=="residue"&Type_detail=="residue"&Diet=="C"&Day=="12"&Replicate=="3"))

  
  # Residue, household waste, control without BSFL - only replictae which had a different feeding rate in comparison to all other containers, not representative

  d <- subset_samples(d, !(Type_detail=="No_larvae"&Diet=="H"))

  
  # Samples with less than 2000 reads and update reads table, only including taxa with > 0 reads
  
  d = prune_samples(sample_sums(d)>=2000, d)
  d <- prune_taxa(taxa_sums(d) > 0, d)    

    # this removed library prep samples for the residue
  
  # update abundance table
  
  d_norm <- transform_sample_counts(d, function(x) x / sum(x))

  
# Sequencing summary ------------------------------------------------------------

  
  # Summary of all samples
    
    d_samples <- subset_samples(d, !Type=="control_DNA" & !Type =="control_lib_prep")
    
    # remove ZOTU without reads
      
    d_samples <- prune_taxa(taxa_sums(d_samples) > 0, d_samples)
  
    # Number of read counts per taxon
    
    taxa_sums(d_samples)
    
    # Number of read counts per sample
    
    sample_sums(d_samples)
    
    # Descriptive statistics on read counts per sample
    
    summary(sample_sums(d_samples)) 
    
    # Summary of phyloseq object
    
    summarize_phyloseq(d_samples)
 
    

  # Summary of DNA extraction controls
  
    d_controls_DNA <- subset_samples(d, Type=="control_DNA")
    
    d_controls_DNA <- prune_taxa(taxa_sums(d_controls_DNA) > 0, d_controls_DNA)
    
    summary(sample_sums(d_controls_DNA)) 
    
    sample_sums(d_controls_DNA)
    
    summarize_phyloseq(d_controls_DNA)
  


# Rarefaction curves ------------------------------------------------------------

rare_d <- rarecurve(t(otu_table(d_samples)), step=50, cex=0.5)

 # plot samples in different color and save



  col <- c("black", "darkred", "forestgreen", "hotpink", "blue")
  lty <- c("solid", "dashed", "dotdash")
  lwd <- c(1, 2)
  pars <- expand.grid(col = col, lty = lty, lwd = lwd, 
                      stringsAsFactors = FALSE)
  
  Nmax <- sapply(rare_d, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(rare_d, max)

  jpeg("output/rare_curves.jpeg", width = 12, height = 12,units="cm",res = 1200)
  

  
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
       ylab = "Species", type = "n")

  for (i in seq_along(rare_d)) {
    N <- attr(rare_d[[i]], "Subsample")
    with(pars, lines(N, rare_d[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
  }

  dev.off()

  # Alpha diversity ---------------------------------------------------------

  # raw reads

  alpha_samples <- plot_richness(d_samples, measures=c("Observed","Shannon","Chao1","Simpson"))


  alpha_samples <-
  
  alpha_samples$data 



  # mean and sd (n>2)
  
  alpha_mean <-
  
  alpha_samples %>% 

    group_by(Diet,Day,Type,Type_detail,variable) %>% 
    
    summarise(n=n(),
              mean=mean(value),
              sd=sd(value))
  
  write.xlsx(alpha_mean , file="output/alpha_mean.xlsx")
  
  
  # diff between samples for n=2
  
  alpha_samples_diff <-
    
    alpha_mean %>% 
    
    filter(n==2) %>% 
    
    select(Diet,Day,Type,Type_detail) %>% 
      
    full_join(alpha_samples,by = c("Diet","Day","Type","Type_detail")) %>% 
    
    group_by(Diet,Day,Type,Type_detail,variable) %>% 
      
    select(-se) %>% 
    
    mutate(diff=lag(value)) %>% 
      
    mutate(diff=value-diff) %>% 
    
    select(Diet,Day,Type,Type_detail,BSFL_detail,value,diff,variable) %>% 
    
    mutate(diff=round(diff,1))
  
  
  write.xlsx(alpha_samples_diff, file="output/alpha_diff.xlsx")

  
# Beta diversity ----------------------------------------------------------


  day <- c("0","3","6","9","12")
   
    
# remove ZOTU without reads from normalized phyloseq ZOTU table
    
    d_norm <- prune_taxa(taxa_sums(d_norm) > 0, d_norm)
          
# samples vs. controls (DNA & library prep)    
  
    
   # only consider taxa that account for at least 1% of relative abundance  
    
    d_filtered <- phyloseq::filter_taxa(d_norm ,function(x) sum(x) > 0.01,TRUE)
      
    
   # NDMS based on UniFrac distance  
    

   d_filtered_ordinate <- ordinate(d_filtered, "NMDS", "unifrac", weighted=TRUE)
    
    
   # Plot ordination  
    
       plot_ordination(d_filtered, d_filtered_ordinate,type="Others",color="Diet",
                    title="") + 
      
      geom_point(size=4,alpha=.5) +
      
      theme_bw(base_size = 15) + 
      
      scale_colour_discrete(name  ="Substrate",
                            breaks=c("C", "Control","H","S"),
                            labels=c("Canteen waste", "Controls","Household waste","Sterile canteen waste"))
      

  # Bacterial compostion of controls is different to samples. But taxon could be possible cross-contaminants. 
   
    
# canteen vs. household
  
  # subset canteen waste and household waste sampels with BSFL
      
  d_norm_C_H <-  subset_samples(d_norm, (Type =="residue" & (Diet=="H"|Diet=="C") &!Type_detail=="No_larvae"))
  
   
  # remove taxa without relative abundance
  
  d_norm_C_H <- prune_taxa(taxa_sums(d_norm_C_H) > 0, d_norm_C_H)
  
  
  # only consider taxa that account for at least 1% of relative abundance  
   
  d_norm_C_H_top <- phyloseq::filter_taxa(d_norm_C_H ,function(x) sum(x) > 0.01,TRUE)
   
  
  # NDMS based on UniFrac distance  
    
  sample_data(d_norm_C_H_top)$Day <- factor(x=sample_data(d_norm_C_H_top)$Day, levels=day)
   
  get_variable(d_norm_C_H_top, "Day")
   
  d_norm_C_H_top_ordinate <- ordinate(d_norm_C_H_top, "NMDS", "unifrac", weighted=TRUE)
   
  
  # Plot ordination  
  
  plot_ordination(d_norm_C_H_top, d_norm_C_H_top_ordinate,type="Others",color="Diet",
                   title="") + 
     
  geom_point(size=4,alpha=.5) +
     
  geom_text(mapping = aes(label = Day), size = 3, vjust = 1.5) +
     
  theme_bw(base_size = 15) + 
        
  scale_colour_manual(values=c("#117733", "#888888"),
                      name  ="Substrate (color)",
                      breaks=c("C", "H"),
                      labels=c("Canteen waste", "Household waste")) +
     
  annotate(geom="text", x=0.2, y=-0.4, label="Stress: 0.06",
                 color="black",fontface=3)

                  
  ggsave("output/beta_C_H.jpeg", height = 12, width = 18, units = 'cm',dpi = "print")
   

#HEEB-treated vs. untreated canteen waste
  
      
  # subset HEEB-treated canteen waste and household waste sampels with BSFL
  
  d_norm_C_S <-  subset_samples(d_norm, Type =="residue" & !Diet=="H")
 
  
  # remove taxa without relative abundance
  
  d_norm_C_S <- prune_taxa(taxa_sums(d_norm_C_S) > 0, d_norm_C_S)
  
   
  # only consider taxa that account for at least 1% of relative abundance  
      
  d_norm_C_S_top  <- phyloseq::filter_taxa(d_norm_C_S,function(x) sum(x) > 0.01,TRUE)
  
 
  # NDMS based on UniFrac distance  
  
  sample_data(d_norm_C_S_top)$Day <- factor(x=sample_data(d_norm_C_S_top)$Day, levels=day)
      
  get_variable(d_norm_C_S_top,"Day")
      
  d_norm_C_S_top_ordinate <- ordinate(d_norm_C_S_top, "NMDS", "unifrac", weighted=TRUE)
  
  # Plot ordination  

  plot_ordination(d_norm_C_S_top, d_norm_C_S_top_ordinate,type="samples",shape="BSFL_detail",color="Diet",
                      title="") + 
        
  geom_point(size=4,alpha=.5) +
        
  geom_text(mapping = aes(label = Day), size = 3, vjust = 1.5) +
        
  theme_bw(base_size = 15) + 
        
  scale_shape_manual(name  ="Sample type (shape)",
                     values=c(17, 15)) +
        
  annotate(geom="text", x=-.32, y=-0.4, label="Stress: 0.05",
                 color="black",fontface=3) +
    
  scale_color_manual(values=c("#117733", "#CC6677"),
                     name  ="Substrate (color)",
                     breaks=c("C", "S"),
                     labels=c("Canteen waste", "Sterile canteen waste"))
      
  ggsave("output/beta_C_S.jpeg", height = 12, width = 18, units = 'cm',dpi = "print")
  

  

#dbRDA ----
  
  # only on canteen and household waste residue samples
  
  #factorize variables
  
  days <- c("0","3","6","9","12")
  replicates <- c("1","2","3","4")
  diets <- c("C","H","S")
  
  
  #prepare data
  
  #convert phyloseq object to tibble
  
  
    #relative abundance
    
  
    ra_tibble <-
      
      as(otu_table(d_norm_C_H_top),"matrix") %>% 
      
      t %>% 
      
      data.frame() %>% 
      
      rownames_to_column(var="sample") 
    
    
    #meta data
  
    
    meta_tibble <-
    
      as(sample_data(d_norm_C_H_top),"matrix") %>% 
    
      data.frame() %>% 
    
      rownames_to_column(var="sample")
  
  
    #merge abundance and meta data
  
    
    ra_meta_tibble <-
      
      left_join(ra_tibble,meta_tibble,by=c("sample")) 
  
  
    #tax data
    
    
    tax_tibble <-
      
      as(tax_table(d_norm_C_H_top),"matrix") %>%
      
      data.frame() %>% 
      
      rownames_to_column(var="ZOTU")
  
  
    #merge abundance, meta, tax data
    
    
    ra_meta_tibble_narrow <- 
      
      ra_meta_tibble %>% 
      
      gather(2:63,key = "ZOTU",value = "count") %>% 
      
      left_join(tax_tibble,by=c("ZOTU")) %>% 
      
      select(sample, Diet, Replicate,Day, Type, Type_detail, ZOTU, count, Phylum, Family,Genus,Species,Class)

        

  #physical-chemical residue compostion (see physical_chemical_temp.R)
  
  
  psychem_rda <-
    
    psychem %>% 
    
    filter(type=="larvae") %>% 
    
    #filter(biowaste=="H"|biowaste=="C") %>% 
    
    rename(Diet = biowaste) %>% 
    
    rename(Day = day) %>% 
    
    rename(Replicate=replicate) %>% 
    
    group_by(Diet,type,parameter) %>% 
    
    spread(key = parameter,value = value) %>% 
    
    ungroup(type) %>% 
    
    select(-type) %>% 
    
    mutate(Day=factor(Day,levels = days)) %>% 
    
    mutate(Replicate=factor(Replicate,levels=replicates)) %>% 
    
    mutate(Diet=factor(Diet,levels = diets))
  
  
  
  #abundance data
  
  
  ra_meta_tibble_narrow_rda <-
    
    ra_meta_tibble_narrow %>% 
    
    select(1:8) %>% 
    
    spread(key = ZOTU,value = count,fill=0) %>% 
    
    as_tibble() %>% 
    
    mutate(Replicate=as.character(Replicate)) %>% 
    
    mutate(Replicate=factor(Replicate,levels=replicates)) %>% 
    
    mutate(Day=case_when(Day=="0"~"0",Day=="3"~"3",Day=="6"~"6",
                         Day=="9"~"9",TRUE~"12")) %>% 
    
    mutate(Diet=factor(Diet,levels = diets)) %>% 
    
    mutate(Day=factor(Day,levels = days))
  
  
  # performance data (see treatment_performance.R)
  
  
  performance_rda <-
    
    performance_messy %>% 
    
    filter(Diet=="C"|Diet=="H") %>% 
    
    filter(!Type=="Control C") %>% 
    
    select(2,4,5,8:12,15) %>% 
    
    left_join(BCR,by=c("Day","Diet","Replicate","DM_Larvae_mg")) %>% 
    
    left_join(larval_protein_accumulation, by=c("Day","Diet","Replicate")) %>% 
    
    select(-N_Larvae_percDM,-C_Larvae_percDM,-Ash_Larvae_percDM,-VS_Larvae_percDM,-Number_larvae_in,
           -Feed_in_gDM) %>% 
    
    mutate(Day=factor(Day,levels = days)) %>% 
    
    mutate(Replicate=factor(Replicate,levels=replicates)) %>% 
    
    mutate(Diet=factor(Diet,levels = diets))
  
  
  # merge the data sets   
  
  
  ra_meta_psychem <- 
    
    ra_meta_tibble_narrow_rda %>% 
    
    left_join(psychem_rda,by=c("Day","Diet","Replicate")) %>% 
    
    left_join(performance_rda,by=c("Day","Diet","Replicate")) %>% 
    
    # correct aw, C_N and N value for day 3
    
    mutate(C_N=case_when(Diet=="C"&Day=="3"~15.2,TRUE~C_N)) %>% 
    
    mutate(N=case_when(Diet=="C"&Day=="3"~3.5,TRUE~N)) %>% 
    
    mutate(C_N=case_when(Diet=="H"&Day=="3"~15.8,TRUE~C_N)) %>% 
    
    mutate(N=case_when(Diet=="H"&Day=="3"~3.3,TRUE~N)) %>%
    
    mutate(aw=round(aw,2)) 
  
  
  # Assess co-lineraity of parameters
  
  performance_psychem <-
  
  performance_C_H %>% 
    
    left_join(psychem_C_H,by=c("Day","Diet","Replicate")) %>% 
    
    select(-1,-2,-3) %>% 
    
    na.omit()
 
   
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(.5, .5, txt, cex = cex.cor * r)
  }
  
  
  jpeg("output/cor.jpeg", width = 30, height = 30,units = "cm",res = 1200)

  pairs(performance_psychem,lower.panel = panel.cor)
 
  dev.off()
  

  # Remove parameters for analysis
  
  ra_meta_psychem <-
  
  ra_meta_psychem %>% 
  
    # ash, VS_ash - correlation with VS
    
    select(-ash,-VS_ash) %>% 
    
    # remove C_N (C:N ratio) - correlation with C, N, pH
    
    select(-C_N) %>% 
  
    # remove waste reduction, bioconversion rate, protein accumulation - correlate with larval weight
    
    select(-Waste_reduction_DM,-bioconversion_rate,-P_gDM_container) %>% 
    
    # remove MC - correlates with volatile solids and C and performance metrics
    
    select(-MC) %>% 
    
    # remove VS - correlates with performance metrics
    
    select(-VS) %>% 
    
    na.omit()
    
 
  #  Is there further multicollinearity between parameters
  
  ev_rda_vif <- 
      
    ra_meta_psychem %>% 
    
    select(aw,N,pH,C,DM_Larvae_mg) %>% 
    
    data.frame() %>% 
    
    na.omit()
  
  
    vif(ev_rda_vif)
    
    vifcor(ev_rda_vif,th = 0.9)
    
    vifstep(ev_rda_vif,th=5) 
  

    # Parameters used in dbRDA without major multicollinearity
    
    # N
    # C
    # aw
    # pH
    # larval weight mg DM
    

  # split data sets again
  
  
  # relative abundance 
  
  
  ra_rda <- 
    
    ra_meta_psychem %>% 
    
    select(7:68) 


  # physico-chemical residue characteristics and larval weight
  
  ev_rda <-   
    
    ra_meta_psychem %>% 
    
    select(1:6,69:73) %>% 
    
    mutate(Diet=as.character(Diet)) %>% 
    
    mutate(Diet=case_when(Diet=="H"~"Household waste",Diet=="C"~"Canteen waste",TRUE~Diet))
  

  
  # normalize & center physical-chemical residue characteristics and larval weight
  
  ev_rda_num <- 
    
    ev_rda %>% 
    
    group_by(Diet,Replicate,Day) %>% 
    
    na.omit() %>% 
    
    gather(7:11,key = parameter,value = value) %>% 
    
    group_by(parameter) %>% 
    
    mutate(value=scale(value,center = TRUE)) %>% 
    
    ungroup() %>% 
    
    spread(key = parameter,value = value) 
    
    
    ev_rda_num_select <-
      
    ev_rda_num %>% 
    
    select(DM_Larvae_mg,N,pH,aw,C)
  
  
  # indentify dissimilary matrix to use 
  
  rankindex(ev_rda_num, ra_rda, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")
  
  
  # run dbRDA analysis  
  
  dbRDA <- capscale(ra_rda ~ aw + C + N + pH + DM_Larvae_mg, ev_rda_num  , method = "bray", add =TRUE)
  
  vector_labels <- c(aw ="aw",C="C",N="N",pH="pH",DM_Larvae_mg="Larval\nweight\nmg DM")
  
  
  # visualize results
  
  ggord(dbRDA,grp_in = ev_rda$Diet, grp_title = "Substrate",addcol="white",vec_lab=vector_labels,ylims=c(-.6,1.2),
        addsize=0.0000000000000,size=ev_rda_num$Day,sizelab="Rearing day",alpha=.7,ellipse = FALSE) + 
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    theme_bw(base_size = 15) 
  

  # check if model is significant
  
  anova(dbRDA)    
  
 
  # check if axis were significant
  
  anova(dbRDA, by="axis", perm.max=1000) 
  
  
  # check environmental parameters that were significant
  
  anova(dbRDA, by="terms", permu=1000)
  
  
  # save plot
  
  ggsave("output/dbRDA.jpeg", height = 12, width = 16, units = 'cm',dpi = "print")
  