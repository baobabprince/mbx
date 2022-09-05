library (tidyverse)
library(dplyr)
library(randomForest)
library(AUC)
library(ggplot2)
setwd('~/Dropbox/Metabolomics')

cv <- function(x){sd(x)/mean(x)}

data_source <- 'MZPos_N1_filtered'   #'_MZPos-N1_biom2'   #'_MZNeg-N1_biom2'    #'MZM10_N2_filtered'
samples_type <- 'Serum'
#samples_type <- 'Stool'
check_val = 'IIRN_Flare'
na_str = c('no_data','_','NA','unknown', 'other','na')

map_file <- paste0(samples_type,'_metadata.txt')
data_file <- paste0(samples_type,'_Data_',data_source,".tsv")  ## For untargeted files
# data_file <- paste0(samples_type,'_Data_',data_source,"_calour.tsv")  ## For targeted files

map0 <- paste0('./Final_datasets/',map_file) %>% read_tsv()
map0$SampleID <- make.names(map0$SampleID)

data1 <- paste0('./Final_datasets/',data_file) %>% read.table(sep="\t", na.strings = na_str, row.names = 1) %>% t() %>% as.data.frame()
#colnames(data1)[1]  <-  'SampleID'
data1 <- data1 %>% rename(SampleID = CompoundID)
data1$SampleID <- make.names(data1$SampleID)
data1 <- data1 %>% mutate_at(c(colnames(data1[-1])), as.numeric) #change the types of the columns (with the peak area values) from character to numeric.

map0 <- map0 %>% filter(SampleID != 'X013_SMC_IIRN_p8s1')  #'X015_SMC_IIRN_p8s3'
map1 <- map0 %>% inner_join(data1 %>% select('SampleID'))

map1 <- map1 %>% filter(Cohort == 'IIRN') %>% filter(Disease_Status == 'CD_Remission')

map1_nf <- map1 %>% filter(IIRN_Flare == 'never_flared')
map1_nf2 <- map1 %>% filter(IIRN_Flare == 'never_flared') %>% filter(Visit == 1)
missing_pn <- setdiff(map1_nf$pn_ID%>% unique(), map1_nf2$pn_ID %>% unique())

map1_nf_misspn <- map1_nf %>% filter(map1_nf$pn_ID %in% missing_pn) %>% filter(Visit == 3)
map1_nf2 <- map1_nf2 %>% rbind(map1_nf_misspn) ## Should be 14 samples(14 subjects)

map1_f <- map1 %>% filter(IIRN_Flare == 'flared')
map1_f2 <- map1 %>% filter(IIRN_Flare == 'flared') %>% filter(Sample_Flare == 'Pre-Flare_3')

missing_pnf <- setdiff(map1_f$pn_ID%>% unique(), map1_f2$pn_ID %>% unique())
map1_f_misspn <- map1_f %>% filter(map1_f$pn_ID %in% missing_pnf) %>% filter(Sample_Flare == 'Pre-Flare_6')
map1_f2 <- map1_f2 %>% rbind(map1_f_misspn) ## Should be 11 samples(11 subjects)
map1_f2$pn_ID %>% n_distinct()

map3 <- map1_nf2 %>% rbind(map1_f2)
map3$pn_ID %>% n_distinct() ## Should be total 25 samples(25 subjects)

map3  <- droplevels(map3)

data2 <- data1 %>% inner_join(map3 %>% select('SampleID'))
#data2$IIRN_Flare <- data2$IIRN_Flare %>% as.factor()
setdiff(map3$SampleID, data2$SampleID)

# pl1 <- 
#   data2 %>% select(-SampleID) %>% apply(2, cv) %>% data.frame() %>% 
#   ggplot(aes(.)) + geom_histogram(color="black", fill = "navy", binwidth = 0.25) + 
#   ggtitle("Serum compounds SD/mean distribution") +
#   #scale_x_continuous(breaks = seq(0, 10, 0.5)) +
#   #scale_y_continuous(breaks = seq(0, 100, 25)) +
#   theme_classic() +
#   #theme(plot.title = element_text(size = 6, hjust = 0.5)) +
#   #theme(axis.text.x = element_text(size = 6)) +
#   theme(aspect.ratio = 0.8)
# pl1
# ggsave(plot = pl1, filename = "Serum_compounds_IIRN1_cv1.png")
# 
data2t <- data2 %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame()
data2t <- data2t %>% mutate(sd1mean = data2t %>% apply(1, cv))
#q25 <- quantile(data$sd1mean)[2]

top75 <- slice_max(data2t, prop=.5, order_by=sd1mean)
##top75 %>% select(-sd1mean) %>% write_delim('Serum_IIRN1_top75compounds.tsv', delim = '\t')

# data3 <- top75 %>% select(-sd1mean) %>% t() %>% as.data.frame()
# data3 <- data3 %>% mutate(SampleID = rownames(data3) ) %>% relocate(SampleID, .before = 1) %>% inner_join(map3 %>% select('SampleID','IIRN_Flare'))

data3 <- data2 %>% inner_join(map3 %>% select('SampleID','IIRN_Flare'))
data3$IIRN_Flare <- data3$IIRN_Flare %>% as.factor()

# pl2 <- 
#   data3 %>% select(-SampleID,-IIRN_Flare) %>% apply(2, cv) %>% data.frame() %>% 
#   ggplot(aes(.)) + geom_histogram(color="black", fill = "navy", binwidth = 0.25) + 
#   ggtitle("Serum compounds SD/mean distribution") +
#   #scale_x_continuous(breaks = seq(0, 10, 0.5)) +
#   #scale_y_continuous(breaks = seq(0, 100, 25)) +
#   theme_classic() +
#   #theme(plot.title = element_text(size = 6, hjust = 0.5)) +
#   #theme(axis.text.x = element_text(size = 6)) +
#   theme(aspect.ratio = 0.8)
# pl2
# ggsave(plot = pl2, filename = "Serum_compounds_IIRN1_cv1_top75.png")

df <- data.frame(1:100, row.names = 1:100)
colnames(df) <- 'auc'
df2 <- data.frame(feature = colnames(data3 %>% select(-SampleID,-IIRN_Flare)))

#set.seed(4)

for (i in 1:100) {
  output.forest <- randomForest(y=data3$IIRN_Flare,
                                ntree  = 1000,
                                x = data3 %>% select(-SampleID, -IIRN_Flare)
                                , na.action = na.omit)
  res = roc(output.forest$votes[,1],factor(1 * (output.forest$y==levels(output.forest$y)[1] )))
  #print( auc(res,min = 0, max = 1) )
  df[i,1] <- auc(res,min = 0, max = 1)
  imp = as.data.frame(importance(output.forest,type = 2))
  imp = data.frame( feature = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini)
  colnames(imp)[-1] <- paste0('MeanDecreaseGini_',as.character(i))
  df2 <- df2 %>% full_join(imp)
}

print(paste0('Mean = ', mean(df$auc)))
print(paste0('SD = ',sd(df$auc)))
#print(paste0('Median = ',median(df$auc)))
#print(paste0('IQR = ',IQR(df$auc)))

# df2t <- df2 %>% t() %>% as.data.frame()
# colnames(df2t) <- df2$feature
# df2t <- df2t[-1,]

df2 <- df2 %>% mutate(impmean = df2 %>% select(-feature) %>% apply(1,mean))
df2 <- df2[rev(order(df2$impmean)), ]

df3 <- df2[1:15,c('feature','impmean')] %>% as.data.frame()

varImpPlot(output.forest,  sort = T, n.var=15, main="Top 15 - Variable Importance",cex=0.6)




res = roc(output.forest$votes[,1],factor(1 * (output.forest$y==levels(output.forest$y)[1] )))
print( auc(res,min = 0, max = 1) )
auc_res = auc(res,min = 0, max = 1)
df = data.frame(pred = output.forest$votes[,1], survived = output.forest$y)

# temp = output.forest$importance
# imp_df = data.frame(asv = row.names(temp), MeanDecreaseGini = as.numeric(temp[,1]))
# ord = rev(order(imp_df$MeanDecreaseGini))
# imp_df = imp_df[ord,]
# 
# top_num = 20
# imp_df_top = imp_df[1:top_num,]
# imp_fig = ggplot(imp_df_top, aes(x=MeanDecreaseGini, y=asv)) + 
#   geom_point() + theme_bw() + ylab('metabolit')
# imp_fig

auc_val = auc(res,min = 0, max = 1)
auc_val = round(auc_val, 3)

res3 = data.frame(fpr = res$fpr, tpr = res$tpr)
p_roc <- ggplot(res3) + 
  geom_line(aes(fpr,tpr)) + xlab("FPR") + ylab("TPR") + theme_bw() + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed") + 
  ggtitle(paste("AUC = ",auc_val)) +
  theme(plot.title = element_text(size = 10))
p_roc
#ggsave(sprintf('%s/%s_auc%f.tiff', out_path, name, auc_res),plot = p_roc, device = 'tiff', width = 3,height =3, compression  = 'lzw')

