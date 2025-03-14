library(tidyverse)
library(data.table)

#work
plates_org <- read.delim("C:/Users/Julius/OneDrive - bwstaff/Promotion/Data/Competition Assay/Species 8_11_23/Counting on Plates 8_11_23.txt")
#home
plates_org <- read.delim("C:/Users/Julius Hoffmann/OneDrive - bwstaff/Promotion/Data/Competition Assay/Species 8_11_23/Counting on Plates 8_11_23.txt")

colnames(plates_org)[c(2,3,5,6,9,11,12,13)] <- c("no_Pf","no_Ec","freq_Pf_plate","freq_Ec_plate","rel_fitness","time.point","stained_species","sample") 

plates_org <- plates_org %>% 
  dplyr::select(-Dilution, -cell.concentration, -doublings, -rel_fitness)

plates_org<- plates_org %>% 
  group_by(sample) %>% 
  arrange(time.point, .by_group = TRUE) %>% 
  mutate(rel_fitness_plate = (log(freq_Pf_plate / lag(freq_Pf_plate, default = first(freq_Pf_plate))) - log(freq_Ec/lag(freq_Ec_plate, default = first(freq_Ec_plate))))/2 ,
         freq_Pf_previous_plate = lag(freq_Pf_plate, default = first(freq_Pf_plate)) ) %>% 
  ungroup()



write_tsv(plates_org,"C:/Users/Julius Hoffmann/OneDrive - bwstaff/Promotion/Data/Competition Assay/Species 8_11_23/Freq_plates_comp_species_8_11_23.tsv" )
