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
info$surv.months <- info$surv.days / 30.4

#parameters and plot
plotdata <- inner_join(info, combine, by = c('patID' = 'bcr'))
sfit <- surv_fit(Surv(surv.months, vitalstatus) ~ ROBO1_level, data = plotdata)
survmedian <- surv_median(sfit)[ ,1:2] # get time for median survival probability
pval <- surv_pvalue(sfit)[2] # pvalue for km estimator
highEGFR_ROBO1table <- cbind.fill(survmedian, pval, fill=NA) # generate output median and pval
write.csv(highEGFR_ROBO1table, file = "highEGFR_ROBO1.csv")
survp <- ggsurvplot(fit = sfit, data = plotdata, 
                    size = 0.8, #line size 
                    palette = c("red", "blue"), #line color: red <- high
                    xlab = 'Time(months)', font.x = 12, break.time.by = 20, xlim = c(0, 120), #x scale
                    ylab = "Survival Probabiliy", font.y = 12, surv.scale = c("percent"), #y scale
                    font.tickslab = c(12, "plain", "plain"), 
                    legend.title = "", legend.labs = c("ROBO1 high", "ROBO1 low"), legend = c(0.85, 0.85), font.legend = 12 #legend
                    # surv.median.line = "hv"
) 
survp <- survp[1]
survp
ggsave(file = "highEGFR-ROBO1.pdf", width = 6, height = 4, 
       units = "in")

