#rpkm table & subset patient barcode and ROBO1 expression level
rpkm <- rnaseq.LUSC[["merged.dat"]]
ROBO1 <- rpkm %>% select(bcr,ROBO1) %>% as.data.frame(., stringsAsFactors = FALSE)
EGFR <- rpkm %>% select(bcr,EGFR) %>% as.data.frame(., stringsAsFactors = FALSE)

#plot ROBO1 expression distribution
library(reshape2)
mROBO1 <- melt(ROBO1, id = 'bcr')
pROBO1 <- ggplot(data = mROBO1[mROBO1$variable == 'ROBO1', ], aes(x = variable, y = value))
pROBO1 + geom_boxplot() + geom_jitter()

#plot EGFR expression distribution
library(reshape2)
mEGFR <- melt(EGFR, id = 'bcr')
pEGFR <- ggplot(data = mEGFR[mEGFR$variable == 'EGFR', ], aes(x = variable, y = value))
pEGFR + geom_boxplot() + geom_jitter(width = 0.2)

#combine robo1 and EGFR expression level
combine <- merge(ROBO1, EGFR)

#expression level of two genes
#robo1
combine$ROBO1_level <- 'na'
combine[combine$ROBO1 >= quantile(combine$ROBO1)[3],'ROBO1_level'] <- 'high'
combine[combine$ROBO1 < quantile(combine$ROBO1)[3],'ROBO1_level'] <- 'low'
#EGFR
combine$EGFR_level <- 'na'
combine[combine$EGFR >= quantile(combine$EGFR)[3],'EGFR_level'] <- 'high'
combine[combine$EGFR < quantile(combine$EGFR)[3],'EGFR_level'] <- 'low'


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
plotdata <- plotdata[plotdata$EGFR_level == 'high', ]
sfit <- surv_fit(Surv(surv.years, vitalstatus) ~ ROBO1_level, data = plotdata)
ggsurvplot(fit = sfit, data = plotdata, risk.table = TRUE)
surv_pvalue(sfit, method = 'FH_p=1_q=1')
surv_pvalue(sfit)


