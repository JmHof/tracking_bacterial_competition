library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(here)

################################
##Evaluation Based on Features##
################################

##Batch Procedure##

path_f_1 <- c("Competition Exp 17_5_23/features/TP0_TP1")
path_f_2 <- c("Competition Exp 17_5_23/features/TP2_TP3")
path_f_3 <- c("Competition Exp 17_5_23/features/TP4_TP5")

plates_samples <- list()

for (p in 1:3) {
all_samples<-data.frame()
all_samples <- (paste("path_f_",p, sep=""))


plates_samples[[p]] <- all_samples  
  
}

path_f_1 <- c("Competition Exp 17_5_23/features/TP0_TP1")
path_f_2 <- c("Competition Exp 17_5_23/features/TP2_TP3")
path_f_3 <- c("Competition Exp 17_5_23/features/TP4_TP5")

plates_samples <- list()


for (p in 1:3) {

#load text files (feature files) from working directory (project directory!)
  
Feature_Files <-
  list.files(path = get(paste("path_f_",p, sep="")),pattern = "*.txt") %>% 
  lapply(., function(x) read.delim2(paste(get(paste("path_f_",p, sep="")),x, sep = "/"), skip = 3 )) 

# get the file names
F_Files <- list.files(path = get(paste("path_f_",p, sep="")),pattern = "*.txt")

#extract certain features (e.g. Intensity_MC_Ch11, Intensity_MC_Ch02) from files and name them with the corresponding sample name (data.frame with dim 3000 x 3).
#Then align this data frame with similar data frames from the other samples.
#Resulting with a data frame containing extracted features for all samples/files
all_samples<-data.frame()
for(t in 1:length(F_Files)){
  
  
  extr_feat <-   Feature_Files[[t]] %>% 
    select(Intensity_MC_Ch11)  #select features to extract
  
  sample <- rep(gsub(".txt", "", F_Files[t]), times = length(extr_feat))  #extract sample name and create vector of same length as extr_feat
  
  
  aligned_features <- cbind(sample,extr_feat) # align sample names with feature columns
  
  #merge data frames from different files/samples
  if(t == 1){
    all_samples<-aligned_features
  }else{
    all_samples <- rbind(all_samples,aligned_features)  
  } 
}

all_samples[,3] <- rep(p,dim(all_samples)[1])

plates_samples[[p]] <- all_samples 

}


features_all_samples <- plates_samples[[1]]

features_all_samples<- base::rbind(features_all_samples,plates_samples[[2]])

features_all_samples<- base::rbind(features_all_samples,plates_samples[[3]])

colnames(features_all_samples) <- c("file_name", "fluorescence_intensity", "plate")

##

##Annotate Data##

features_all_samples_an <- features_all_samples %>% 
  separate_wider_position(file_name ,widths = c(row= 1,  well =2), cols_remove = FALSE)%>%
  as.data.frame()  %>% 
  mutate(well = as.numeric(well)) %>% 
  mutate(time = case_when(well <=6 & plate == 1 ~ 0,
                                well >=7 & plate == 1 ~ 2,
                                well <=6 & plate == 2 ~ 4,
                                well >=7 & plate == 2 ~ 6,
                                well <=6 & plate == 3 ~ 8,
                                well >=7 & plate == 3 ~ 10),
         sample = case_when(well == 1 ~ 1,
                            well == 2 ~ 2,
                            well == 3 ~ 3,
                            well == 4 ~ 4,
                            well == 5 ~ 5,
                            well == 6 ~ 6,
                            well == 7 ~ 1,
                            well == 8 ~ 2,
                            well == 9 ~ 3,
                            well == 10 ~ 4,
                            well == 11 ~ 5,
                            well == 12 ~ 6)) %>% 
  mutate(stained_isolate = case_when(row == "A" & sample == 1 ~ "I6_stained_ctrl",
                                     row == "A" & sample == 2 ~ "I6_unstained",
                                     row == "A" & sample == 3 ~ "blank",
                                     row == "A" & sample == 4 ~ "pred_ctrl",
                                     row == "A" & sample == 5 ~ "I7_stained_ctrl",
                                     row == "A" & sample == 6 ~ "I7_unstained",
                                     row == "B" ~ "I7",
                                     row == "C" ~ "I7",
                                     row == "D" ~ "I7",
                                     row == "E" ~ "I6",
                                     row == "F" ~ "I6",
                                     row == "G" ~ "I6",
                                     row == "H" & sample == 1 ~ "I6_stained_ctrl",
                                     row == "H" & sample == 2 ~ "I6_unstained",
                                     row == "H" & sample == 3 ~ "blank",
                                     row == "H" & sample == 4 ~ "pred_ctrl",
                                     row == "H" & sample == 5 ~ "I7_stained_ctrl",
                                     row == "H" & sample == 6 ~ "I7_unstained"),
         predator = case_when(
           row == "A" & sample == 4 ~ "YES",
           row == "A" & sample != 4 ~ "NO",
           row == "B" & sample <= 3 ~ "YES",
           row == "B" & sample >= 4 ~ "NO",
           row == "C" & sample <= 3 ~ "YES",
           row == "C" & sample >= 4 ~ "NO",
           row == "D" & sample <= 3 ~ "YES",
           row == "D" & sample >= 4 ~ "NO",
           row == "E" & sample <= 3 ~ "YES",
           row == "E" & sample >= 4 ~ "NO",
           row == "F" & sample <= 3 ~ "YES",
           row == "F" & sample >= 4 ~ "NO",
           row == "G" & sample <= 3 ~ "YES",
           row == "G" & sample >= 4 ~ "NO",
           row == "H" & sample == 4 ~ "YES",
           row == "H" & sample != 4 ~ "NO"
         ),
         starting_freq_I7 = case_when(
           row == "A" & sample %in% 1:4 ~ 0,
           row == "A" & sample %in% 5:6 ~ 1,
           row == "H" & sample %in% 1:4 ~ 0,
           row == "H" & sample %in% 5:6 ~ 1,
           
           !row %in% c("A","H") & sample %in% c(1,4) ~ 0.9,
           !row %in% c("A","H") & sample %in% c(2,5) ~ 0.5,
           !row %in% c("A","H") & sample %in% c(3,6) ~ 0.1
         )
  )

##calculate mean fluorescence in cultures

features_all_samples_an <- features_all_samples_an %>% 
  group_by(row, well, time, sample, stained_isolate, predator, starting_freq_I7) %>% 
  reframe(mean_intensity = mean(fluorescence_intensity)) %>% 
  mutate(culture_ID = paste(row,sample,sep = "_"), .before= 1) 

##Visualising mean fluorescence##

starting.labs <- c("t(0) 10vol% Pf","t(0) 50vol% Pf","t(0) 90vol% Pf" )
names(starting.labs) <- c(0.1,0.5, 0.9)

stain.labs <- c("I6 stained", "I7 stained", "I7 unstained", "I7 stained control", "I6 unstained", "I6 stained control")
names(stain.labs) <- c("I6", "I7", "I7_unstained", "I7_stained_ctrl", "I6_unstained", "I6_stained_ctrl")

##Plot fluorescence Intensity over time for monocultures

features_all_samples_an %>% 
  dplyr::filter(starting_freq_I7 == 1 | starting_freq_I7 == 0) %>% 
  dplyr::filter(!stained_isolate == c("blank"), !stained_isolate == c("pred_ctrl")) %>% 
  ggplot( aes(x = time, y = mean_intensity, colour = as.factor(stained_isolate), group = culture_ID))+
  geom_line(linewidth=0.3)+
  geom_point(size=0.5)+
  facet_grid( ~ stained_isolate , labeller = labeller( stained_isolate = stain.labs, starting_freq_I7 = starting.labs)  )+
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

ggsave("I6_I7_monoculture_fluorescence.jpg", device = "jpg", path = "Results/", width = 350, height = 150,  units = "mm" )

##Visualising Retention Curves##

#before retitions curves can be visualised the fraction of stained cells per sample needs to be added to the data. 
#Although this could be newly calculated based on the features, I will instead use the already calculated values from another data file.

stain_freqs_isolate <- read.delim(here("Competition Exp 17_5_23", "Freq_ISX_comp_17_5_23__th_19_2_25.tsv"))
stain_freqs_isolate <-stain_freqs_isolate %>% 
mutate(time = time.point*2)

features_all_samples_an <- left_join(features_all_samples_an, stain_freqs_isolate, by = join_by(row,well,time,starting_freq_I7, predator, sample),suffix = c("", ".y"))

#Plot retention curves for stained mono- and co-cultures. For these curves mean fluorescence of samples is divided by the fraction of stained cells to provide a comparable measuore of fluorescence intensity between mono- and co-cultures. 

features_all_samples_an %>% 
  dplyr::filter(!stained_isolate == c("blank"),!stained_isolate == c("I6_unstained"), !stained_isolate == c("I7_unstained"), !stained_isolate == c("pred_ctrl")) %>% 
  ggplot( aes(x = time, y = mean_intensity/ratio, colour = as.factor(starting_freq_I7), group = culture_ID))+
  geom_line(linewidth=0.3)+
  geom_point(size=0.5)+
  facet_grid( ~ stained_isolate , labeller = labeller( stained_isolate = stain.labs, starting_freq_I7 = starting.labs)  )+
  scale_color_discrete(name = "starting frequency I7")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  #ylim(0,1000)+
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

ggsave("I6_and_I7_Fluorescence_retention_curves.jpg", device = "jpg", path = "Results/", width = 350, height = 150,  units = "mm" )
