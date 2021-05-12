library(tidyverse)
library(magrittr)
data <- read_csv("Serum1_Data_TF.csv")
data <- read_csv("TF-N1_206s.csv")

data %<>% select(-contains("149"), -contains("150"), -contains("151"))
metadata <- read_csv("Serum1_Metadata.csv")
groups <- read_csv("metabolites groups.csv")

names(groups)[2] <- "Compound"

all <- 
  data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  left_join(groups) %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(Compound %in% groups$Compound) 

pl <- 
  ggplot(all, aes(x = Disease_Status, y = value, fill = Disease_Status)) + 
  geom_boxplot() + 
  # geom_jitter(alpha = .3)+
  facet_wrap(~Compound, scales = "free_y") +
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()
        , strip.text = element_text(size = 6), strip.background = element_blank())

ggsave(paste0("all_with_fasting", "_metabolites.png"), pl, height = 20, width = 30)
