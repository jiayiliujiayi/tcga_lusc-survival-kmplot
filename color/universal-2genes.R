setwd(dir = "/Users/JiayiLiu/Downloads/ngsanalysis/color")
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
DRIV <- rpkm %>% select(bcr,DRIV) %>% as.data.frame(., stringsAsFactors = FALSE)
SECOND <- rpkm %>% select(bcr,SECOND) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine DRIV and SECOND expression level
combine <- merge(SECOND, DRIV)

#plot DRIV expression distribution
mDRIV <- melt(DRIV, id = 'bcr')
pDRIV <- ggplot(data = mDRIV[mDRIV$variable == 'DRIV', ], aes(x = variable, y = value))
pDRIV + geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  labs(title = 'overall DRIV expression level distribution', 
       x = '', 
       y = 'Expression Level')

#DRIV high or low
DRIVq <- quantile(combine$DRIV) %>% as.vector()
combine[combine$DRIV >= DRIVq[3], 'DRIV_level'] <- 'high'
combine[combine$DRIV < DRIVq[3], 'DRIV_level'] <- 'low'
combine <- combine[combine$DRIV_level == 'high', ]

#plot SECOND expression distribution
mSECOND <- melt(combine %>% select(bcr, SECOND), id = 'bcr')
pSECOND <- ggplot(data = mSECOND[mSECOND$variable == 'SECOND', ], aes(x = variable, y = value))
pSECOND + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high DRIV, SECOND expression level distribution', 
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
info$surv.months <- info$surv.days / 30.4

#ggplot
plotdata <- inner_join(info, combine, by = c('patID' = 'bcr'))
sfit <- surv_fit(Surv(surv.months, vitalstatus) ~ SECOND_level, data = plotdata)
surv_pvalue(sfit)
ggsurvplot(fit = sfit, data = plotdata, 
           palette = c("red", "blue"), 
           risk.table = TRUE, 
           title = 'high DRIV, SECOND KM plot; DRIV cutoff: 50%, SECOND cutoff: 50%, pval = 0.039'
) 



