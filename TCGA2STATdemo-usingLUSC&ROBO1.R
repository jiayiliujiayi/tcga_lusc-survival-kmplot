setwd(dir = "/Users/JiayiLiu/Downloads/ngsanalysis")
# using TCGA2STAT package
library(TCGA2STAT)
library(dplyr)
library(survival)
library(survminer)

#import rpkm table and clincal information
rnaseq.LUSC <- getTCGA(disease="LUSC", 
                       data.type="RNASeq", 
                       type="RPKM", 
                       clinical = TRUE )
#rpkm table & subset patient barcode and robo1 expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
ROBO1 <- rpkm %>% select(bcr,ROBO1) %>% as.data.frame(., stringsAsFactors = FALSE)
SLIT2 <- rpkm %>% select(bcr,SLIT2) %>% as.data.frame(., stringsAsFactors = FALSE)

#plot robo1 expression distribution
library(reshape2)
mROBO1 <- melt(ROBO1, id = 'bcr')
p1 <- ggplot(data = mROBO1[mROBO1$variable == 'ROBO1', ], aes(x = variable, y = value))
p1 + geom_boxplot() + geom_jitter()

# quantile(ROBO1$ROBO1)
# 0%        25%        50%        75%       100% 
# 0.4583946  3.2356763  6.0379079  9.5345882 74.8973593 
ROBO1$ROBO1_level <- 'na'
ROBO1[ROBO1$ROBO1 >= 6.0379079,'ROBO1_level'] <- 'high'
ROBO1[ROBO1$ROBO1 < 6.0379079,'ROBO1_level'] <- 'low'

#clinical information & add patient ID column
clinical <- rnaseq.LUSC[['clinical']] %>% as.data.frame(., stringsAsFactors = FALSE)
clinical$patID <- rownames(clinical)


#subset clinical information
info <- clinical %>% select(patID, vitalstatus, daystodeath, daystolastfollowup)
info$surv.days <- ifelse(info$vitalstatus == 1, info$daystodeath, info$daystolastfollowup)
info$vitalstatus <- info$vitalstatus %>% as.numeric()
info$surv.days <- info$surv.days %>% as.numeric()

#ggplot
plotdata <- inner_join(info, ROBO1, by = c('patID' = 'bcr'))
sfit <- surv_fit(Surv(surv.days, vitalstatus) ~ ROBO1_level, data = plotdata)
ggsurvplot(fit = sfit, data = plotdata, risk.table = TRUE)
surv_pvalue(sfit, method = 'FH_p=1_q=1')

# variable       pval   method  pval.txt
# ROBO1_level 0.06074044 Log-rank p = 0.061
