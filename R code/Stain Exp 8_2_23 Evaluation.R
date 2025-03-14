library(tidyverse)
library(data.table)
library(patchwork)
library(grid)


################################
##Evaluation Based on Features##
################################

##Batch Procedure##

#load text files (feature files) from working directory (project directory!)
path_f <- c("Stain Concentration Experiment 8_2_23/Features")
Feature_Files <-
  list.files(path = "Stain Concentration Experiment 8_2_23/Features",pattern = "*.txt") %>% 
  lapply(., function(x) read.delim2(paste(path_f,x, sep = "/"), skip = 3 )) 

# get the file names
F_Files <- list.files(path = "Stain Concentration Experiment 8_2_23/Features",pattern = "*.txt")

#extract certain features (e.g. Intensity_MC_Ch11, Intensity_MC_Ch02) from files and name them with the corresponding sample name (data.frame with dim 3000 x 3).
#Then align this data frame with similar data frames from the other samples.
#Resulting with a data frame containing extracted features for all samples/files
all_samples<-data.frame()
for(t in 1:length(F_Files)){
  
  
extr_feat <-   Feature_Files[[t]] %>% 
    select(Intensity_MC_Ch11,Intensity_MC_Ch02)  #select features to extract

sample <- rep(gsub(".txt", "", F_Files[t]), times = length(extr_feat))  #extract sample name and create vector of same length as extr_feat
 

aligned_features <- cbind(sample,extr_feat) # align sample names with feature columns

#merge data frames from different files/samples
  if(t == 1){
          all_samples<-aligned_features
          }else{
          all_samples <- rbind(all_samples,aligned_features)  
  } 
 
}



all_sample_annotated <- all_samples %>% 
   separate_wider_position(sample ,widths = c(row= 1,  well =2)) %>% 
   mutate(well = as.numeric(well)) %>% 
   mutate(stained_isolate = case_when( well == 1  ~ "I8",
                                       well == 3 ~ "I8",
                                       well == 5 ~ "I8",
                                       well == 7 ~ "I8",
                                       well == 9 ~ "I8",
                                       well == 11 ~ "I8",
                                       well == 2 ~ "I7",
                                       well == 4 ~ "I7",
                                       well == 6 ~ "I7",
                                       well == 8 ~ "I7",
                                       well == 10 ~ "I7",
                                       well == 12 ~ "I7"),
          time_point= case_when( well == 1 | well == 2 ~ 0,
                                 well == 3 | well == 4 ~ 1,
                                 well == 5 | well == 6 ~ 2,
                                 well == 7 | well == 8 ~ 3,
                                 well == 9 | well == 10 ~ 4,
                                 well == 11 | well == 12 ~ 5),
          stain_concentration = case_when(row == "B" ~ "Ctrl_1",
                                          row == "C" ~ "400µ",
                                          row == "D" ~ "300µ",
                                          row == "E" ~ "200µ",
                                          row == "F" ~ "100µ",
                                          row == "G" ~ "Ctrl_2",))


str(all_sample_annotated)   
all_sample_annotated$stained_isolate <- as.factor(all_sample_annotated$stained_isolate)
all_sample_annotated$stain_concentration <- as.factor(all_sample_annotated$stain_concentration)
all_sample_annotated$row <- as.factor(all_sample_annotated$row)
all_sample_annotated$well <- as.factor(all_sample_annotated$well)

write_tsv(all_sample_annotated, "Stain Concentration Experiment 8_2_23/Features_extracted.tsv" )

##Ch11 Fluorescence over time

all_sample_annotated %>% 
  group_by(stained_isolate, stain_concentration, time_point) %>% 
  reframe(mean_intensity = mean(Intensity_MC_Ch11)) %>% 

ggplot()+
  geom_line(aes(x = time_point,y = mean_intensity, colour = interaction(stained_isolate, stain_concentration)))+
  scale_colour_discrete(name = "stained isolate + \nstain volume")+
  #scale_y_continuous(trans= "log2")+
  facet_grid(~stained_isolate)+
  labs(title = "Ch11 intensity TP0-5")+
  theme_bw()+
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 20, face = "bold"),
         legend.text = element_text(size = 18),
         legend.title = element_text(size = 18, face = "bold"),
         legend.key.size = unit(1, "cm"),
         axis.text = element_text(color = "grey50", size = 18),
         axis.title = element_text(size = 20),
         axis.title.y = element_text( colour = "black"),
         axis.title.y.right = element_text( colour = "#00A9E0"),
         strip.text = element_text(size = 18, colour = "white", face = "bold"),
         strip.background = element_rect(color="black", fill="#008ECE", size=1, linetype="solid"))

##Boxplots of Ch11 fluorescence over time 
ggplot(data = all_sample_annotated, aes(x = interaction(stained_isolate, stain_concentration) , y=Intensity_MC_Ch11+500 ))+
  geom_boxplot(aes(fill = interaction(stained_isolate, stain_concentration)))+
  #scale_fill_manual(name = "sample",values=c("blue","red","orange","pink","yellow","green"), labels=c("Ctrl1","400","300","200","100","Ctrl2"))+
  scale_y_continuous(trans= "log2")+
  labs(title = "Ch11 intensity at TP0")+
  facet_grid(~stained_isolate)
  theme_bw()+
  theme( plot.title = element_text(size = 20, face = "bold"),
         legend.text = element_text(size = 18),
         legend.title = element_text(size = 18, face = "bold"),
         legend.key.size = unit(1, "cm"),
         axis.text = element_text(color = "grey50", size = 18),
         axis.title = element_text(size = 20),
         axis.title.y = element_text( colour = "#FF8E7B"),
         axis.title.y.right = element_text( colour = "#00A9E0"),
         strip.text = element_text(size = 18, colour = "white", face = "bold"),
         strip.background = element_rect(color="black", fill="#008ECE", size=1, linetype="solid"))




##Bscatterplot from samples at the same time point and isolates
ggplot(data = all_sample_annotated, aes(x = Intensity_MC_Ch11+400 , y=Intensity_MC_Ch02 ))+
  geom_point(aes(colour = interaction(stained_isolate, stain_concentration)))+
  scale_colour_discrete(name = "sample")+
  scale_y_continuous(trans= "log2", limits = c(1,35000))+
  scale_x_continuous(trans= "log2",limits = c(1,300000))+
  labs(title = "Ch11 vs Ch2 intensity")+
  facet_grid(stained_isolate~time_point)
  theme_bw()+
  theme( plot.title = element_text(size = 20, face = "bold"),
         legend.text = element_text(size = 18),
         legend.title = element_text(size = 18, face = "bold"),
         legend.key.size = unit(1, "cm"),
         axis.text = element_text(color = "grey50", size = 18),
         axis.title = element_text(size = 20),
         axis.title.y = element_text( colour = "#FF8E7B"),
         axis.title.y.right = element_text( colour = "#00A9E0"),
         strip.text = element_text(size = 18, colour = "white", face = "bold"),
         strip.background = element_rect(color="black", fill="#008ECE", size=1, linetype="solid"))

 
##save plots
#home
ggsave("Ch11 vs Ch2 intensity after 0h I8.jpeg", device = "jpeg", path = "Stain Concentration Experiment 8_2_23/Evaluation",  dpi = 300, width = 17, height = 11,  units = "in" )  


####################################
##Evaluation Based on Thresholding##
####################################


stain_freqs <- read.csv("Stain Concentration Experiment 8_2_23/Stain Concentration Experiment 8_2_23_stain_counts_thresholding.csv")


ISX_freqs <- stain_freqs[1:84,] %>% 
  separate_wider_delim(file, delim = "." , names = c("sample", "drop"), cols_remove = TRUE) %>% 
  select(sample ,ratio, stained, total) %>% 
  separate_wider_position(sample ,widths = c(bs= 1,  sample =3)) %>% 
  select(sample ,ratio, stained, total) %>% 
  separate_wider_position(sample ,widths = c(row = 1,  well =2)) %>% 
  dplyr::filter(!row == "A" & !row == "H" ) %>% 
  mutate(well = as.numeric(well))
  

ISX_freqs_annotated <- ISX_freqs %>% 
  mutate(stained_isolate = case_when( well == 1  ~ "I8",
                                      well == 3 ~ "I8",
                                      well == 5 ~ "I8",
                                      well == 7 ~ "I8",
                                      well == 9 ~ "I8",
                                      well == 11 ~ "I8",
                                      well == 2 ~ "I7",
                                      well == 4 ~ "I7",
                                      well == 6 ~ "I7",
                                      well == 8 ~ "I7",
                                      well == 10 ~ "I7",
                                      well == 12 ~ "I7"),
         time_point= case_when( well == 1 | well == 2 ~ 0,
                                well == 3 | well == 4 ~ 1,
                                well == 5 | well == 6 ~ 2,
                                well == 7 | well == 8 ~ 3,
                                well == 9 | well == 10 ~ 4,
                                well == 11 | well == 12 ~ 5),
         stain_concentration = case_when(row == "B" ~ "Ctrl_1",
                                         row == "C" ~ "400µ",
                                         row == "D" ~ "300µ",
                                         row == "E" ~ "200µ",
                                         row == "F" ~ "100µ",
                                         row == "G" ~ "Ctrl_2",))
ISX_freqs_annotated %>% 

ggplot(aes(x = time_point, y = ratio))+
  geom_line(aes( colour = interaction(stained_isolate, stain_concentration)))+
  facet_grid(~stained_isolate)+
  scale_colour_discrete(name = "stained isolate + \nstain volume")+
  labs(y= "ratio stained/unstained")+
  ylim(0,1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

