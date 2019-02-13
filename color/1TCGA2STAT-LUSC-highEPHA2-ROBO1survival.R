setwd(dir = "/Users/JiayiLiu/Downloads/ngsanalysis/color")
library(TCGA2STAT)
library(dplyr)
library(survival)
library(survminer)
library(reshape2)
library(TCGA2STAT)
library(rowr)

# import rpkm table and clincal information
rnaseq.LUSC <- getTCGA(disease="LUSC", 
                       data.type="RNASeq", 
                       type="RPKM", 
                       clinical = TRUE )
load("rnaseq-LUSC.RData")
#rpkm table & subset patient barcode and expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
EPHA2 <- rpkm %>% select(bcr,EPHA2) %>% as.data.frame(., stringsAsFactors = FALSE)
ROBO1 <- rpkm %>% select(bcr,ROBO1) %>% as.data.frame(., stringsAsFactors = FALSE)
#combine EPHA2 and ROBO1 expression level
combine <- merge(ROBO1, EPHA2)

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

#plot ROBO1 expression distribution
mROBO1 <- melt(combine %>% select(bcr, ROBO1), id = 'bcr')
pROBO1 <- ggplot(data = mROBO1[mROBO1$variable == 'ROBO1', ], aes(x = variable, y = value))
pROBO1 + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  labs(title = 'high EPHA2, ROBO1 expression level distribution', 
       x = '', 
       y = 'Expression Level')

#expression level of two genes
#ROBO1
ROBO1q <- quantile(combine$ROBO1) %>% as.vector()
combine[combine$ROBO1 >= ROBO1q[3], 'ROBO1_level'] <- 'high'
combine[combine$ROBO1 < ROBO1q[3], 'ROBO1_level'] <- 'low'

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

#parameters and plot
plotdata <- inner_join(info, combine, by = c('patID' = 'bcr'))
sfit <- surv_fit(Surv(surv.months, vitalstatus) ~ ROBO1_level, data = plotdata)
survmedian <- surv_median(sfit)[ ,1:2] # get time for median survival probability
pval <- surv_pvalue(sfit)[2] # pvalue for km estimator
highEPHA2_ROBO1table <- cbind.fill(survmedian, pval, fill=NA) # generate output median and pval
write.csv(highEPHA2_ROBO1table, file = "highEPHA2_ROBO1.csv")
survp <- ggsurvplot(fit = sfit, data = plotdata, 
           size = 0.8, #line size 
           palette = c("red", "blue"), #line color: red <- high
           xlab = 'Time(months)', font.x = 12, break.time.by = 20, xlim = c(0, 120), #x scale
           ylab = "Survival", font.y = 12, surv.scale = c("percent"), #y scale
           font.tickslab = c(12, "plain", "plain"), legend = "none"
           #legend.title = "", legend.labs = c("ROBO1 high", "ROBO1 low"), legend = c(0.85, 0.85), font.legend = 12 #legend
           # surv.median.line = "hv"
           ) 
survp <- survp[1]
survp
 ggsave(file = "highEPHA2-ROBO1.pdf", width = 6, height = 4, 
       units = "in")

