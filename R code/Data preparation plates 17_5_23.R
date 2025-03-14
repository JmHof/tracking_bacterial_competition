library(tidyverse)
library(data.table)

plates_org <- read.delim("Competition Exp 17_5_23/Counting on plates 17_5_23.txt")


plates_org <- plates_org[,1:11] %>% 
  mutate(cell_concentration = case_when( time.point == 0 ~ (1/(0.1*(1/5)^(dilution-1)*0.1)*total)*(800/1300),
                                          time.point != 0 ~ (1/(0.1*(1/5)^(dilution-1)*0.1)*total))) %>% 
  mutate(density_I6 = freq_I6 * cell_concentration,
         density_I7 = freq_I7 * cell_concentration) %>% 
  group_by(replicate) %>% 
  arrange(time.point, .by_group = TRUE) %>% 
  mutate(rel_fitness = (log(density_I7 / lag(density_I7, default = first(density_I7))) - log(density_I6/lag(density_I6, default = first(density_I6))))/2 ,
         freq_I7_previous = lag(freq_I7, default = first(freq_I7)) ) %>% 
  ungroup()



write_tsv(plates_org,"Competition Exp 17_5_23/Freq_on_plates_comp_17_5_23_analysis.tsv" )
