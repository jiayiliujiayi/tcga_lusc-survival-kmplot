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
EFNA1 <- rpkm %>% select(bcr,EFNA1) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine EPHA2 and EFNA1 expression level
combine <- merge(EFNA1, EPHA2)

#plot EPHA2 expression distribution
mEPHA2 <- melt(EPHA2, id = 'bcr')
pEPHA2 <- ggplot(data = mEPHA2[mEPHA2$variable == 'EPHA2', ], aes(x = variable, y = value))
pEPHA2 + geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  labs(title = 'overall EPHA2 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#EPHA2 high or low
EPHA2q <- quantile(combine$EPHA2) %>% as.vector()
combine[combine$EPHA2 >= EPHA2q[3], 'EPHA2_level'] <- 'high'
combine[combine$EPHA2 < EPHA2q[3], 'EPHA2_level'] <- 'low'
combine <- combine[combine$EPHA2_level == 'high', ]

#plot EFNA1 expression distribution
mEFNA1 <- melt(combine %>% select(bcr, EFNA1), id = 'bcr')
pEFNA1 <- ggplot(data = mEFNA1[mEFNA1$variable == 'EFNA1', ], aes(x = variable, y = value))
pEFNA1 + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high EPHA2, EFNA1 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#EFNA1
EFNA1q <- quantile(combine$EFNA1) %>% as.vector()
combine[combine$EFNA1 >= EFNA1q[3], 'EFNA1_level'] <- 'high'
combine[combine$EFNA1 < EFNA1q[3], 'EFNA1_level'] <- 'low'

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
sfit <- surv_fit(Surv(surv.years, vitalstatus) ~ EFNA1_level, data = plotdata)
surv_pvalue(sfit)
ggsurvplot(fit = sfit, data = plotdata, 
           risk.table = TRUE, 
           title = 'high EPHA2, EFNA1 KM plot; EPHA2 cutoff: 50%, EFNA1 cutoff: 50%, pval = 0.095'
) 



