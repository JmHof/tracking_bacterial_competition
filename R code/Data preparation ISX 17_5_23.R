library(tidyverse)
library(data.table)
library(ggplot2)
library(rstatix)
library(lmerTest)
#library(lme4)


##Processing export data from ISX/Ole for 17_5_23 Assay##
#########################################################
#stains_org <- read.csv("Competition Exp 17_5_23/Competition Assay 17_5_23_ISX_thresholding_org.csv")

##Ole Johannsen determined thresholds for ioslates I6 and I7 separately using ROCs. He then applied these thresholds on the Fluorescence intensitiy values in Ch11 extracted from the ISX data. 
##He did this by applying each of the two thresholds on the complete data set. This resulted in two files each one for each threshold.
##Thus, in each file half of the samples was thresholded with the wrong threshold (e.g. the I7 threshold although I6 was stained). 
##Therefore, these files need to be combined by extracting those samples for which the appropriate threshold was used. This mean that if isolate I6 was stained then this sample needs to be extracted from the file in which the I6 threshold was applied.
##This is done by first annotating the files and then selecting the correct rows of both data sets. Additionally, controls where bacteria where present are separately extracted from the file with the I7 threshold.

stains_I7 <- read.csv("Competition Exp 17_5_23/Ole/Thresholds_19_2_25/stain_counts_thresholding_Intensity_MC_Ch11_201.35.csv")
stains_I6 <- read.csv("Competition Exp 17_5_23/Ole/Thresholds_19_2_25/stain_counts_thresholding_Intensity_MC_Ch11_858.49.csv")


##for loop to annotate both data sets
for (i in c(1:2)) {

isolate <- c("I7","I6")  
stains_org <- get(paste("stains_",isolate[i], sep=""))  

 
ISX_freqs <- stains_org %>% 
      separate_wider_regex(file, c(time.points = ".*", "/", well = ".*")) %>% 
      separate_wider_delim(well, delim = "." , names = c("well", "drop"), cols_remove = TRUE) %>% 
      select(time.points,well,ratio, stained, total) %>% 
      separate_wider_position(well ,widths = c(row= 1,  well =2)) %>% 
      separate_wider_delim(time.points, delim = " & " , names = c("first_TP", "second_TP"), cols_remove = TRUE) %>% 
      separate_wider_position(second_TP ,widths = c(drop= 2,  TP =2)) %>% 
      select(TP,well,row, ratio, stained, total)

ISX_freqs$well <- as.numeric(ISX_freqs$well)  
ISX_freqs$TP <- as.numeric(ISX_freqs$TP) 


ISX_freqs <- ISX_freqs %>% 
  mutate(time.point = case_when(well <=6 & TP == 1 ~ 0,
                                well >=7 & TP == 1 ~ 1,
                                well <=6 & TP == 3 ~ 2,
                                well >=7 & TP == 3 ~ 3,
                                well <=6 & TP == 5 ~ 4,
                                well >=7 & TP == 5 ~ 5),
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
   mutate(stained_isolate = case_when(row == "A" & sample == 1 ~ 6,
                                     row == "A" & sample == 2 ~ 6,
                                     row == "A" & sample == 3 ~ 0,
                                     row == "A" & sample == 4 ~ 0,
                                     row == "A" & sample == 5 ~ 7,
                                     row == "A" & sample == 6 ~ 7,
                                     row == "B" ~ 7,
                                     row == "C" ~ 7,
                                     row == "D" ~ 7,
                                     row == "E" ~ 6,
                                     row == "F" ~ 6,
                                     row == "G" ~ 6,
                                     row == "H" & sample == 1 ~ 6,
                                     row == "H" & sample == 2 ~ 6,
                                     row == "H" & sample == 3 ~ 0,
                                     row == "H" & sample == 4 ~ 0,
                                     row == "H" & sample == 5 ~ 7,
                                     row == "H" & sample == 6 ~ 7),
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
            row == "A" & sample ==  1 ~ 0,
            row == "A" & sample ==  2 ~ 0,
            row == "A" & sample ==  3 ~ 0,
            row == "A" & sample ==  4 ~ 0,
            row == "A" & sample ==  5 ~ 1,
            row == "A" & sample ==  6 ~ 1,
            row == "H" & sample ==  1 ~ 0,
            row == "H" & sample ==  2 ~ 0,
            row == "H" & sample ==  3 ~ 0,
            row == "H" & sample ==  4 ~ 0,
            row == "H" & sample ==  5 ~ 1,
            row == "H" & sample ==  6 ~ 1,
            row != c("A","H") & sample == c(1) ~ 0.9,
            row != c("A","H") & sample == c(4) ~ 0.9,
            row != c("A","H") & sample == c(2) ~ 0.5,
            row != c("A","H") & sample == c(5) ~ 0.5,
            row != c("A","H") & sample == c(3) ~ 0.1,
            row != c("A","H") & sample == c(6) ~ 0.1
          )) %>% 
  select(-TP)

ISX_freqs$stained_isolate <- as.numeric(ISX_freqs$stained_isolate) 
str(ISX_freqs)

#create "freq_I7" from "ratio"
ISX_freqs[,"freq_I7"]<- ISX_freqs[,"ratio"]

#select rows of samples where I6 was stained and calculate the correct I7 frequency
ISX_freqs[grep(6,ISX_freqs$stained_isolate), "freq_I7"] <-   (1-filter(ISX_freqs,stained_isolate== 6)[,"ratio"])

#create freq_I6 and calculate I6 frequency from I7 frequency
ISX_freqs <- ISX_freqs %>% 
  mutate(freq_I6 = 1- freq_I7)

ISX_freqs <- ISX_freqs %>% 
  group_by(row,sample) %>% 
  arrange(time.point, .by_group = TRUE) %>% 
  mutate(rel_fitness = (log(freq_I7 / lag(freq_I7, default = first(freq_I7))) - log(freq_I6/lag(freq_I6, default = first(freq_I6))))/2 ,
         freq_I7_previous = lag(freq_I7, default = first(freq_I7)),
         replicate = paste(row, sample, sep = "")) %>% 
  ungroup()

#each annotated data set is saved
write_tsv(ISX_freqs, paste("Competition Exp 17_5_23/Ole/Thresholds_19_2_25/threshold_",isolate[i],".tsv", sep=""))

}
#write_tsv(ISX_freqs, "Competition Exp 17_5_23/Freq_ISX_comp_17_5_23_analysis.tsv" )

##The annotated files are loaded, appropriate data is extracted and then combined in a single file.
##This file is then used for the analysis

ISX_freqs_I7 <- read.delim("Competition Exp 17_5_23/Ole/Thresholds_19_2_25/threshold_I7.tsv")
ISX_freqs_I6 <- read.delim("Competition Exp 17_5_23/Ole/Thresholds_19_2_25/threshold_I6.tsv")

no_stained <- ISX_freqs_I7 %>% 
                      filter(stained_isolate == 0)#these are controls without bacteria i.e. blank, only predators
I7_stained <- ISX_freqs_I7 %>% 
                      filter(stained_isolate == 7)
I6_stained <- ISX_freqs_I6 %>% 
                      filter(stained_isolate == 6)

ISX_freqs_combined <- rbind(no_stained,I7_stained,I6_stained)

write_tsv(ISX_freqs_combined, "Competition Exp 17_5_23/Freq_ISX_comp_17_5_23__th_19_2_25.tsv" )

##correction for fading staining## 
##################################
#work
stains <- read.delim("C:/Users/Julius/OneDrive - bwstaff/Promotion/Data/Competition Assay/17_5_23/Competition Assay 17_5_23 ISX_incl_densities.tsv")

#home
stains <- read.delim("C:/Users/Julius Hoffmann/OneDrive - bwstaff/Promotion/Data/Competition Assay/17_5_23/Competition Assay 17_5_23 ISX_incl_densities.tsv")


String_Patterns = c('A', 'H')

stain_corr <- stains %>% 
  filter(str_detect(replicate, str_c(String_Patterns, collapse = "|"))) %>% 
  select(replicate, total, stained, stained_isolate, predator, time.point, freq_I7, freq_I6, density_I7, density_I6,Objects_mL)

String_Patterns_2 = c('A1', 'H1','A7', 'H7')
String_Patterns_3 = c('A5', 'H5','A10', 'H10')

stain_fit_I6 <- stain_corr %>% 
  filter(str_detect(replicate, str_c(String_Patterns_2, collapse = "|"))) 

stain_fit_I7 <- stain_corr %>% 
  filter(str_detect(replicate, str_c(String_Patterns_3, collapse = "|"))) 

model <- nls(stained ~ -a * exp(time.point) + c, start = list(a = 10, c = 2000), data = stain_fit_I7)

#assess model fit
sse <- as.vector((summary(model)[[3]])^2*10) #calculate sse by the residual standard error of the model and DF
null <- lm(stained~1, data = stain_fit_I7) #a null model
sst <- as.vector(unlist(summary.aov(null)[[1]][2])) #work out the total variation present in y by comparing to a null model
sst
100*(sst-sse)/sst 

#predict fitness from model
newdata  = data.frame(time.point = c(0,1,2,3,4,5))
pred_stained <- predict(model,newdata)/2000
pred_stained_inv <- 1/pred_stained # correction factors for TP 0-5 and I7
stain_corr_I7 <- data.frame(corr_factor = c(pred_stained_inv,1,1,1,1,1,1) , time.point = c(0,1,2,3,4,5,0,1,2,3,4,5), stained_isolate = c(7,7,7,7,7,7,6,6,6,6,6,6)) #creates a data.frame with correction factors at the time point for I7. 
#It also contains values for I6 so that it can easily be merges with the original data. Thus, far the correction factor for I6 is a dummy value (1).
stain_corr_I7$stained_isolate <- as.factor(stain_corr_I7$stained_isolate)

stain_corrected<- stains %>% 
  filter(!replicate %like% "A") %>% 
  filter(!replicate %like% "H") %>% 
  select(stained, total, replicate, starting_freq_I7, stained_isolate, predator, time.point, Objects_mL)

stain_corrected$stained_isolate <- as.factor(stain_corrected$stained_isolate)

stain_corrected <- right_join(stain_corrected,stain_corr_I7, by =c("time.point", "stained_isolate"), unmatched = "drop" )

stain_corrected <- stain_corrected %>% 
  mutate(stained_corrected = stained * corr_factor) %>% 
  mutate(freq_I7 = stained_corrected/total ,
         ratio = stained_corrected/total)

#select rows of samples where I6 was stained and calculate the correct I7 frequency
stain_corrected[grep(6,stain_corrected[,"stained_isolate"]), "freq_I7"] <-   (1-filter(stain_corrected,stained_isolate== 6)[,"ratio"])

#create freq_I6 and calculate I6 frequency from I7 frequency
stain_corrected <- stain_corrected %>% 
  mutate(freq_I6 = 1- freq_I7)

#calculate cell densities of I6 and I7 stain_corrected <- 
stain_corrected<- stain_corrected %>% 
  mutate(density_I6 = freq_I6 * Objects_mL,
         density_I7 = freq_I7 * Objects_mL) %>% 
  group_by(replicate) %>% 
  arrange(time.point, .by_group = TRUE) %>% 
  mutate(rel_fitness = (log(density_I7 / lag(density_I7, default = first(density_I7))) - log(density_I6/lag(density_I6, default = first(density_I6))))/2 ,
         freq_I7_previous = lag(freq_I7, default = first(freq_I7)) ) %>% 
  ungroup()


write_tsv(stain_corrected, "C:/Users/Julius/OneDrive - bwstaff/Promotion/Data/Competition Assay/17_5_23/Freq_ISX_corr_comp_17_5_23_analysis.tsv" )
