library(tidyverse)
library(data.table)
library(ggplot2)
library(lmerTest)
#library(lme4)
library(rstatix)

##set working directory to user directory
#setwd("..")

##Processing export data from IS for 8_11_23 Assay##
####################################################
#1st run
#plate1_org <- read.delim("OneDrive - bwstaff/Promotion/Data/Competition Assay/Competition Exp Species 8_11_23/Competition_CH11_900_th_TP0-2.txt", skip = 3, dec = ",")
#plate2_org <- read.delim("OneDrive - bwstaff/Promotion/Data/Competition Assay/Competition Exp Species 8_11_23/Competition_CH11_144_th_TP3-5.txt", skip = 3, dec = ",")

#2nd run
#plate1_org <- read.delim("Competition Exp Species 8_11_23//Competition_CH11_420_th_TP0-2_2nd_run.txt", skip = 3, dec = ",")
#plate2_org <- read.delim("Competition Exp Species 8_11_23//Competition_CH11_420_th_TP3-5_2nd_run.txt", skip = 3, dec = ",")

#2nd run with stain dependent thresholds (based on both runs)
#plate1_org <- read.delim("Competition Exp Species 8_11_23/Competition_Ch11_420_680_th_TP0-2_2nd_run.txt", skip = 3, dec = ",")
#plate2_org <- read.delim("Competition Exp Species 8_11_23/Competition_Ch11_420_680_th_TP3-5_2nd_run.txt", skip = 3, dec = ",")

#2nd run with stain dependent thresholds (based on 2nd run)
plate1_org <- read.delim("Promotion/Data/Competition Assay/Competition Exp Species 8_11_23/Competition_Ch11_320_379_th_TP0-2_2nd_run.txt", skip = 3, dec = ",")
plate2_org <- read.delim("Promotion/Data/Competition Assay/Competition Exp Species 8_11_23/Competition_Ch11_320_379_th_TP3-5_2nd_run.txt", skip = 3, dec = ",")


#1st run data plate preparation 
#Plate 1
ISX_freqs_p1 <- plate1_org[1:(dim(plate1_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(well, X.Gated..unstained..All, X.Gated..stained..All) %>% 
  separate_wider_position(well ,widths = c(row= 1,  well =2)) 
ISX_freqs_p1$well <- as.numeric(ISX_freqs_p1$well)
#Plate 2
ISX_freqs_p2 <- plate2_org[1:(dim(plate2_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(well,X...unstained..All,X...stained..All) %>% 
  separate_wider_position(well ,widths = c(row= 1,  well =2)) 
ISX_freqs_p2$well <- as.numeric(ISX_freqs_p2$well)



#2nd run data plate preparation 
#Plate 1
ISX_freqs_p1 <- plate1_org[1:(dim(plate1_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(well, freq_unstained, freq_stained) %>% 
  separate_wider_position(well ,widths = c(row= 1,  well =2)) 
ISX_freqs_p1$well <- as.numeric(ISX_freqs_p1$well)
#Plate 2 
ISX_freqs_p2 <- plate2_org[1:(dim(plate2_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(well,freq_unstained, freq_stained) %>% 
  separate_wider_position(well ,widths = c(row= 1,  well =2)) 
ISX_freqs_p2$well <- as.numeric(ISX_freqs_p2$well)


#Plate 1 sample assignment
ISX_freqs_p1_t <- ISX_freqs_p1 %>% 
  mutate(time.point = case_when(well <=4 ~ 0,
                                well >=5 & well <=8 ~ 1,
                                well >=9  ~ 2,
                                ),
         sample = case_when(well == 1 ~ 1,
                            well == 2 ~ 2,
                            well == 3 ~ 3,
                            well == 4 ~ 4,
                            well == 5 ~ 1,
                            well == 6 ~ 2,
                            well == 7 ~ 3,
                            well == 8 ~ 4,
                            well == 9 ~ 1,
                            well == 10 ~ 2,
                            well == 11 ~ 3,
                            well == 12 ~ 4)) %>% 
  mutate(stained_species = case_when(row == "A" & sample == 1 ~ "Ec",
                                     row == "A" & sample == 2 ~ "Ec",
                                     row == "A" & sample == 3 ~ "none",
                                     row == "A" & sample == 4 ~ "none",
                                  
                                     row == "B" ~ "Pf",
                                     row == "C" ~ "Pf",
                                     row == "D" ~ "Pf",
                                     row == "E" ~ "Ec",
                                     row == "F" ~ "Ec",
                                     row == "G" ~ "Ec",
                                     row == "H" & sample == 1 ~ "Pf",
                                     row == "H" & sample == 2 ~ "Pf",
                                     row == "H" & sample == 3 ~ "none",
                                     row == "H" & sample == 4 ~ "none",
                                     )
         ) %>% 
  mutate(sample = paste(row,sample, sep = ""))

#Plate 2 sample assignment
ISX_freqs_p2_t <- ISX_freqs_p2 %>% 
  mutate(time.point = case_when(well <=4 ~ 3,
                                well >=5 & well <=8 ~ 4,
                                well >=9  ~ 5,
  ),
  sample = case_when(well == 1 ~ 1,
                     well == 2 ~ 2,
                     well == 3 ~ 3,
                     well == 4 ~ 4,
                     well == 5 ~ 1,
                     well == 6 ~ 2,
                     well == 7 ~ 3,
                     well == 8 ~ 4,
                     well == 9 ~ 1,
                     well == 10 ~ 2,
                     well == 11 ~ 3,
                     well == 12 ~ 4)) %>% 
  mutate(stained_species = case_when(row == "A" & sample == 1 ~ "Ec",
                                     row == "A" & sample == 2 ~ "Ec",
                                     row == "A" & sample == 3 ~ "none",
                                     row == "A" & sample == 4 ~ "none",
                                     
                                     row == "B" ~ "Pf",
                                     row == "C" ~ "Pf",
                                     row == "D" ~ "Pf",
                                     row == "E" ~ "Ec",
                                     row == "F" ~ "Ec",
                                     row == "G" ~ "Ec",
                                     row == "H" & sample == 1 ~ "Pf",
                                     row == "H" & sample == 2 ~ "Pf",
                                     row == "H" & sample == 3 ~ "none",
                                     row == "H" & sample == 4 ~ "none",
  )
  ) %>% 
  mutate(sample = paste(row,sample, sep = ""))

colnames(ISX_freqs_p1_t)[3:4] <- c("unstained","stained")
colnames(ISX_freqs_p2_t)[3:4] <- c("unstained","stained")

ISX_freqs <- rbind(ISX_freqs_p1_t, ISX_freqs_p2_t)
  


ISX_freqs$stained <- as.numeric(ISX_freqs$stained) 




#create "freq_I7" from "ratio"
ISX_freqs[,"freq_Pf"]<- ISX_freqs[,"stained"]

#select rows of samples where I6 was stained and calculate the correct I7 frequency
ISX_freqs[grep("Ec",ISX_freqs$stained_species), "freq_Pf"] <-   (100-filter(ISX_freqs, stained_species== "Ec")[,"stained"])

#create freq_I6 and calculate I6 frequency from I7 frequency
ISX_freqs <- ISX_freqs %>% 
  mutate(freq_Ec = 100 - freq_Pf)

#calculate relative fitness
ISX_freqs <- ISX_freqs %>% 
  group_by(row,sample) %>% 
  arrange(time.point, .by_group = TRUE) %>% 
  mutate(rel_fitness_ISX = (log(freq_Pf / lag(freq_Pf, default = first(freq_Pf))) - log(freq_Ec/lag(freq_Ec, default = first(freq_Ec))))/2 ,
         freq_Pf_previous_ISX = lag(freq_Pf, default = first(freq_Pf)) ) %>% 
  ungroup()

colnames(ISX_freqs)[c(8,9)] <- c("freq_Pf_ISX","freq_Ec_ISX")

ISX_freqs <- ISX_freqs %>% 
  mutate(starting_frequency = case_when((well == 1 | well == 2| well == 3| well == 5| well == 6| well == 7| well == 9 | well == 10| well == 11) & row == "A" ~ 0,
                                        (well == 1 | well == 2| well == 3| well == 5| well == 6| well == 7| well == 9 | well == 10| well == 11) & row == "H" ~ 1,
                                        well == 1 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.9,
                                        well == 2 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.5,
                                        well == 3 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.1,
                                        well == 4 & (row == "B" | row == "C") ~ 1,
                                        well == 4 & (row == "F" | row == "G") ~ 0,
                                        well == 5 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.9,
                                        well == 6 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.5,
                                        well == 7 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.1,
                                        well == 8 & (row == "B" | row == "C") ~ 1,
                                        well == 8 & (row == "F" | row == "G") ~ 0,
                                        well == 9 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.9,
                                        well == 10 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.5,
                                        well == 11 & (row == "B" | row == "C"| row == "D"| row == "E"| row == "F"| row == "G") ~ 0.1,
                                        well == 12 & (row == "B" | row == "C") ~ 1,
                                        well == 12 & (row == "F" | row == "G") ~ 0),
         freq_Pf_ISX = freq_Pf_ISX * 0.01,
         freq_Ec_ISX = freq_Pf_ISX * 0.01)

write_tsv(ISX_freqs, "Promotion/Data/Competition Assay/Competition Exp Species 8_11_23/Freq_ISX_comp_species_8_11_23_Ch11_TH_320_379.tsv" )


#plot frequencies derived from stained samples
ISX_freqs %>% 
  filter(!row == "A" & !row == "H" & !well == c(4,8,12)) %>% 
ggplot( aes(x = time.point, y = freq_Pf_ISX))+
  geom_line(aes(linetype = as.factor(stained_species), colour = sample))+
  facet_grid(~ starting_frequency )#+
  ylim(0,1)
xlim(1,5)

ggsave("ISX_frequencies_all_TP_threshold_420_680_Ch11.jpeg", device = "jpeg", path = "C:/Users/Julius/OneDrive - bwstaff/Promotion/Data/Competition Assay/Competition Exp Species 8_11_23",  dpi = 250, width = 17, height = 11,  units = "in" )               
