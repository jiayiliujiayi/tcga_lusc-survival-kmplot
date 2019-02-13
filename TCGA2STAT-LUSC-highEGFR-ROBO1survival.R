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
  #import mutation data
  #mut.LUSC <- getTCGA(disease = "LUSC", 
                    Edata.type = "Mutation", 
                    type = "somatic")
  #mutation table & subset EGFR mutation data
  #mut <- mut.LUSC[["dat"]]
  #Emut <- t(mut) %>% as.data.frame(., stringsAsFactors = FALSE)
  #EGFRmutation <- mut %>% select(EGFR)
  #EGFRmutation$patID <- rownames(EGFRmutation)
  #subset EGFR mutated patients
  #EGFRmut <- EGFRmutation %>% subset(., EGFR == "1")
  
#rpkm table & subset patient barcode and expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
EGFR <- rpkm %>% select(bcr,EGFR) %>% as.data.frame(., stringsAsFactors = FALSE)
ROBO1 <- rpkm %>% select(bcr,ROBO1) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine EGFR and ROBO1 expression level
combine <- merge(ROBO1, EGFR)

#plot EGFR expression distribution
mEGFR <- melt(EGFR, id = 'bcr')
pEGFR <- ggplot(data = mEGFR[mEGFR$variable == 'EGFR', ], aes(x = variable, y = value))
pEGFR + geom_boxplot() + 
  geom_jitter(width = 0.2) + 
  labs(title = 'overall EGFR expression level distribution', 
       x = '', 
       y = 'Expression Level')

#EGFR high or low
EGFRq <- quantile(combine$EGFR) %>% as.vector()
combine[combine$EGFR >= EGFRq[3], 'EGFR_level'] <- 'high'
combine[combine$EGFR < EGFRq[3], 'EGFR_level'] <- 'low'
combine <- combine[combine$EGFR_level == 'high', ]

#plot ROBO1 expression distribution
mROBO1 <- melt(combine %>% select(bcr, ROBO1), id = 'bcr')
pROBO1 <- ggplot(data = mROBO1[mROBO1$variable == 'ROBO1', ], aes(x = variable, y = value))
pROBO1 + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high EGFR, ROBO1 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#ROBO1
ROBO1q <- quantile(combine$ROBO1) %>% as.vector()
combine[combine$ROBO1 >= ROBO1q[4], 'ROBO1_level'] <- 'high'
combine[combine$ROBO1 < ROBO1q[4], 'ROBO1_level'] <- 'low'

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
sfit <- surv_fit(Surv(surv.years, vitalstatus) ~ ROBO1_level, data = plotdata)
surv_pvalue(sfit, method = "GB" )
ggsurvplot(fit = sfit, data = plotdata, 
           risk.table = TRUE, 
           title = 'high EGFR, ROBO1 KM plot; EGFR cutoff: 50%, ROBO1 cutoff: 50%, \n
           pval = 0.17 (method Gehan-Breslow)'
) 
