setwd(dir = "/Users/JiayiLiu/Downloads/ngsanalysis")
# using TCGA2STAT package
library(TCGA2STAT)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

#import rpkm table and clincal information
rnaseq.LUSC <- getTCGA(disease="LUSC", 
                       data.type="RNASeq", 
                       type="RPKM", 
                       clinical = TRUE )
#rpkm table & subset patient barcode and SLIT2 expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
ROBO1 <- rpkm %>% select(bcr,ROBO1) %>% as.data.frame(., stringsAsFactors = FALSE)
SLIT2 <- rpkm %>% select(bcr,SLIT2) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine ROBO1 and SLIT2 expression level
combine <- merge(SLIT2, ROBO1)

#plot ROBO1 expression distribution
library(reshape2)
mROBO1 <- melt(ROBO1, id = 'bcr')
pROBO1 <- ggplot(data = mROBO1[mROBO1$variable == 'ROBO1', ], aes(x = variable, y = value))
pROBO1 + geom_boxplot() + geom_jitter(width = 0.2)  + 
  labs(title = 'overall ROBO1 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#ROBO1 high or low
ROBO1q <- quantile(combine$ROBO1) %>% as.vector()
combine[combine$ROBO1 >= ROBO1q[4], 'ROBO1_level'] <- 'high'
combine[combine$ROBO1 < ROBO1q[4], 'ROBO1_level'] <- 'low'
combine <- combine[combine$ROBO1_level == 'high', ]

#plot SLIT2 expression distribution
library(reshape2)
mSLIT2 <- melt(combine %>% select(bcr, SLIT2), id = 'bcr')
pSLIT2 <- ggplot(data = mSLIT2[mSLIT2$variable == 'SLIT2', ], aes(x = variable, y = value))
pSLIT2 + geom_boxplot() + geom_jitter(width = 0.2)  +
  labs(title = 'high ROBO1, SLIT2 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#SLIT2
SLIT2q <- quantile(combine$SLIT2) %>% as.vector()
combine[combine$SLIT2 >= SLIT2q[2], 'SLIT2_level'] <- 'high'
combine[combine$SLIT2 < SLIT2q[2], 'SLIT2_level'] <- 'low'

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
ggsurvplot(fit = sfit, data = plotdata, 
           risk.table = TRUE, 
           title = 'high ROBO1, SLIT2 KM plot; EPHA2 cutoff: 50%, ROBO1 cutoff: 50%, pval = 0.065'
) 
surv_pvalue(sfit, method = 'FH_p=1_q=1')
surv_pvalue(sfit)
