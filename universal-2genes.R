setwd(dir = "/Users/JiayiLiu/Downloads/ngsanalysis")
library(TCGA2STAT)
library(dplyr)
library(survival)
library(survminer)
library(reshape2)

#import rpkm table and clincal information
rnaseq.LUSC <- getTCGA(disease="LUSC", 
                       data.type="RNASeq", 
                       type="RPKM", 
                       clinical = TRUE )
#rpkm table & subset patient barcode and expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
DRIVING <- rpkm %>% select(bcr,DRIVING) %>% as.data.frame(., stringsAsFactors = FALSE)
SECOND <- rpkm %>% select(bcr,SECOND) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine DRIVING and SECOND expression level
combine <- merge(SECOND, DRIVING)

#plot DRIVING expression distribution
mDRIVING <- melt(DRIVING, id = 'bcr')
pDRIVING <- ggplot(data = mDRIVING[mDRIVING$variable == 'DRIVING', ], aes(x = variable, y = value))
pDRIVING + geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  labs(title = 'overall DRIVING expression level distribution', 
       x = '', 
       y = 'Expression Level')

#DRIVING high or low
DRIVINGq <- quantile(combine$DRIVING) %>% as.vector()
combine[combine$DRIVING >= DRIVINGq[3], 'DRIVING_level'] <- 'high'
combine[combine$DRIVING < DRIVINGq[3], 'DRIVING_level'] <- 'low'
combine <- combine[combine$DRIVING_level == 'high', ]

#plot SECOND expression distribution
mSECOND <- melt(combine %>% select(bcr, SECOND), id = 'bcr')
pSECOND <- ggplot(data = mSECOND[mSECOND$variable == 'SECOND', ], aes(x = variable, y = value))
pSECOND + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high DRIVING, SECOND expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#SECOND
SECONDq <- quantile(combine$SECOND) %>% as.vector()
combine[combine$SECOND >= SECONDq[3], 'SECOND_level'] <- 'high'
combine[combine$SECOND < SECONDq[3], 'SECOND_level'] <- 'low'

#clinical information & add patient ID column
clinical <- rnaseq.LUSC[['clinical']] %>% as.data.frame(., stringsAsFactors = FALSE)
clinical$patID <- rownames(clinical)

#subset clinical information
info <- clinical %>% select(patID, vitalstatus, daystodeath, daystolastfollowup)
info$surv.days <- ifelse(info$vitalstatus == 1, info$daystodeath, info$daystolastfollowup)
info$vitalstatus <- info$vitalstatus %>% as.numeric()
info$surv.days <- info$surv.days %>% as.numeric() 
info$surv.years <- info$surv.days / 365

#ggplot
plotdata <- inner_join(info, combine, by = c('patID' = 'bcr'))
sfit <- surv_fit(Surv(surv.years, vitalstatus) ~ SECOND_level, data = plotdata)
surv_pvalue(sfit)
ggsurvplot(fit = sfit, data = plotdata, 
           risk.table = TRUE, 
           title = 'high DRIVING, SECOND KM plot; DRIVING cutoff: 50%, SECOND cutoff: 50%, pval = 0.095'
) 



