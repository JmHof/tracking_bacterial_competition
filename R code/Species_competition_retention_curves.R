library(tidyverse)
library(data.table)
library(ggplot2)
library(lmerTest)
#library(lme4)
library(rstatix)
library(here)
##set working directory to user directory
#setwd("..")

##Processing export data from IS for 8_11_23 Assay##
####################################################

#Load ISX data
plate1_org <- read.delim(here("Competition Exp Species 8_11_23", "Competition_Ch11_320_379_th_TP0-2_2nd_run.txt"), skip = 3, dec = ",")
plate2_org <- read.delim(here("Competition Exp Species 8_11_23", "Competition_Ch11_320_379_th_TP3-5_2nd_run.txt"), skip = 3, dec = ",")

#2nd run data plate preparation 
#Plate 1
ISX_freqs_p1 <- plate1_org[1:(dim(plate1_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(-drop) %>% 
  separate_wider_position(well ,widths = c(row= 1,  well =2)) 
ISX_freqs_p1$well <- as.numeric(ISX_freqs_p1$well)
#Plate 2 
ISX_freqs_p2 <- plate2_org[1:(dim(plate2_org)[1]-1),] %>% 
  separate_wider_delim(File, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
  select(-drop) %>% 
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

colnames(ISX_freqs)[c(9,10)] <- c("freq_Pf_ISX","freq_Ec_ISX")

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

#correctly annotate stained and unstained monocultures.

ISX_freqs <- ISX_freqs %>% 
mutate(stained_species = case_when(sample == "A1" ~ "Ec_stained_ctrl",
                                   sample == "A2" ~ "Ec_unstained",
                                   sample == "A3" ~ "blank",
                                   
                                   
                                   sample == "B4" ~ "Pf_stained_ctrl",
                                   
                                   
                                   sample == "C4" ~ "Pf_unstained",
                                   
                                   
                                   sample == "F4"~ "Ec_stained_ctrl",
                                   
                                   
                                   sample == "G4"~ "Ec_unstained",
                                   
                                   sample == "H1" ~ "Pf_stained_ctrl",
                                   sample == "H2" ~ "Pf_unstained",
                                   sample == "H3" ~ "blank",
                                   .default = as.character(stained_species)
                                   
))

##Visualising Retention Curves##

starting.labs <- c("t(0) 10vol% Pf","t(0) 50vol% Pf","t(0) 90vol% Pf" )
names(starting.labs) <- c(0.1,0.5, 0.9)

stain.labs <- c("Ec stained", "Pf stained", "Pf unstained", "Pf stained control", "Ec unstained", "Ec stained control")
names(stain.labs) <- c("Ec", "Pf", "Pf_unstained", "Pf_stained_ctrl", "Ec_unstained", "Ec_stained_ctrl")

##Plot fluorescence Intensity over time for monocultures

ISX_freqs %>% 
  dplyr::filter(starting_frequency == 1 | starting_frequency == 0) %>% 
  dplyr::filter(!stained_species == c("blank")) %>% 
  mutate(time = time.point*2) %>%
ggplot( aes(x = time, y = IntensityCh11, colour = as.factor(stained_species), group = sample))+
  geom_line(linewidth=0.3)+
  geom_point(size=0.5)+
  facet_grid( ~ stained_species , labeller = labeller( stained_species = stain.labs, starting_frequency = starting.labs)  )+
  scale_color_discrete(guide = "none")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  labs(x= "time [h]", y= "fluorescence intensity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1, "cm"),
        axis.text = element_text(color = "grey50", size = 18),
        axis.title = element_text(size = 20),
        axis.title.y = element_text( colour = "black"),
        axis.title.y.right = element_text( colour = "#00A9E0"),
        strip.text = element_text(size = 18, colour = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="white", linewidth=1, linetype="solid"))

ggsave("Ec_Pf_monoculture_fluorescence.jpg", device = "jpg", path = "Results/", width = 350, height = 150,  units = "mm" )

#Plot retention curves for stained mono- and co-cultures. For these curves mean fluorescence of samples is divided by the fraction of stained cells to provide a comparable measuore of fluorescence intensity between mono- and co-cultures. 

ISX_freqs %>% 
  dplyr::filter(!stained_species == c("blank"),!stained_species == c("Ec_unstained"), !stained_species == c("Pf_unstained")) %>% 
  mutate(time = time.point*2) %>%
  ggplot( aes(x = time, y = IntensityCh11/stained, colour = as.factor(starting_frequency), group = sample))+
  geom_line(linewidth=0.3)+
  geom_point(size=0.5)+
  facet_grid( ~ stained_species , labeller = labeller( stained_species = stain.labs, starting_frequency = starting.labs)  )+
  scale_color_discrete(name = "starting frequency Pf")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  ylim(0,1000)+
  labs(x= "time [h]", y= "mean fluorescence intensity/\nfraction of stained cells")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1, "cm"),
        axis.text = element_text(color = "grey50", size = 18),
        axis.title = element_text(size = 20),
        axis.title.y = element_text( colour = "black"),
        axis.title.y.right = element_text( colour = "#00A9E0"),
        strip.text = element_text(size = 18, colour = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="white", linewidth=1, linetype="solid"))

ggsave("Ec_and_Pf_Fluorescence_retention_curves.jpg", device = "jpg", path = "Results/", width = 350, height = 150,  units = "mm" )
