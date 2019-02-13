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
EPHA2 <- rpkm %>% select(bcr,EPHA2) %>% as.data.frame(., stringsAsFactors = FALSE)
SLIT2 <- rpkm %>% select(bcr,SLIT2) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine epha2 and SLIT2 expression level
combine <- merge(SLIT2, EPHA2)

#plot epha2 expression distribution
mEPHA2 <- melt(EPHA2, id = 'bcr')
pEPHA2 <- ggplot(data = mEPHA2[mEPHA2$variable == 'EPHA2', ], aes(x = variable, y = value))
pEPHA2 + geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  labs(title = 'overall EPHA2 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#epha2 high or low
EPHA2q <- quantile(combine$EPHA2) %>% as.vector()
combine[combine$EPHA2 >= EPHA2q[3], 'EPHA2_level'] <- 'high'
combine[combine$EPHA2 < EPHA2q[3], 'EPHA2_level'] <- 'low'
combine <- combine[combine$EPHA2_level == 'high', ]

#plot SLIT2 expression distribution
mSLIT2 <- melt(combine %>% select(bcr, SLIT2), id = 'bcr')
pSLIT2 <- ggplot(data = mSLIT2[mSLIT2$variable == 'SLIT2', ], aes(x = variable, y = value))
pSLIT2 + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high EPHA2, SLIT2 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#SLIT2
SLIT2q <- quantile(combine$SLIT2) %>% as.vector()
combine[combine$SLIT2 >= SLIT2q[3], 'SLIT2_level'] <- 'high'
combine[combine$SLIT2 < SLIT2q[3], 'SLIT2_level'] <- 'low'

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
sfit <- surv_fit(Surv(surv.years, vitalstatus) ~ SLIT2_level, data = plotdata)
surv_pvalue(sfit)
ggsurvplot(fit = sfit, data = plotdata, 
           risk.table = TRUE, 
           title = 'high EPHA2, SLIT2 KM plot; EPHA2 cutoff: 50%, SLIT2 cutoff: 50%, pval = 0.11'
) 



