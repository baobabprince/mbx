library(tidyverse)
library(magrittr)
data <- read_csv("Serum1_Data_TF.csv")
data <- read_csv("TF-N1_206s.csv")

data %<>% select(-contains("149"), -contains("150"), -contains("151"))
metadata <- read_csv("Serum1_Metadata.csv")
groups <- read_csv("metabolites groups.csv")

names(groups)[2] <- "Compound"

for(element in unique(groups$group)){
merged_table <- 
data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  left_join(groups) %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(group == element) %>% filter(Fasting == "Y")
}  
pl <- 
 ggplot(merged_table, aes(x = Disease_Status, y = value, fill = Disease_Status)) + 
  geom_boxplot() + 
  # geom_jitter(alpha = .3)+
  facet_wrap(~Compound, scales = "free_y") +
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()
        , strip.text = element_text(size = 6), strip.background = element_blank())
pl
ggsave(paste0(unique(merged_table$group), "_metabolites.png"), pl)
}
###################################################
for(element in unique(groups$group)){
merged_table <- 
  data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  group_by(Compound) %>% 
  summarise(freq = value/sum(value), SampleID, Compound) %>% 
  left_join(groups, by = "Compound") %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(group == element) 

pl <- 
  ggplot(merged_table, aes(x = Disease_Status, y = freq, fill = Disease_Status)) + 
  geom_boxplot() + 
  # geom_jitter(alpha = .3)+
  facet_wrap(~Compound, scales = "free_y") +
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()
        , strip.text = element_text(size = 6), strip.background = element_blank())
ggsave(paste0(unique(merged_table$group), "_freq_metabolites.png"), pl)
}


merged_table <- 
  data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  group_by(Compound) %>% 
  summarise(freq = value/sum(value), SampleID, Compound) %>% 
  left_join(groups, by = "Compound") %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(group == "TCA") 

pl <- 
  ggplot(merged_table, aes(x = Disease_Status, y = freq, fill = Disease_Status)) + 
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(alpha = .3)+
  scale_y_continuous(limits = quantile(c(0, 0.1))) +
  facet_grid(~Compound, scales = "free_y") 

# plotly::ggplotly(pl)


merged_table <- 
  data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  left_join(groups) %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(Compound == "Guanine") %>% filter(Fasting == "Y") %>% 
  filter(value < 3e6)

pl <- 
  ggplot(merged_table, aes(x = Disease_Status, y = value, fill = Disease_Status)) + 
  geom_boxplot() + 
  # geom_jitter(alpha = .3)+
  facet_wrap(~Compound, scales = "free_y") +
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()
        , strip.text = element_text(size = 6), strip.background = element_blank())
pl
ggsave(paste0("Guanine", "_metabolites.png"), pl)


all <- 
data %>% 
  pivot_longer(contains("_"), names_to = "SampleID", values_to = "value") %>% 
  left_join(groups) %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(Fasting == "Y")


pl <- 
  ggplot(all, aes(x = Disease_Status, y = value, fill = Disease_Status)) + 
  geom_boxplot() + 
  # geom_jitter(alpha = .3)+
  facet_wrap(~Compound, scales = "free_y") +
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()
        , strip.text = element_text(size = 6), strip.background = element_blank())
pl
ggsave(paste0("all", "_metabolites.png"), pl, height = 30, width = 40, units = "cm")
