setwd("/Users/JiayiLiu/Downloads/ngsanalysis/")

library(dplyr);library(stringr); library(ggplot2); 
#import ROBO1 expression level, add patient barcode
originalROBO1 <- read.delim("ROBO1.tsv", stringsAsFactors=FALSE)
originalROBO1$site <- 'na'
originalROBO1$site <- substr(originalROBO1$sample, start = 14, stop = 15)

ROBO1 <- originalROBO1[originalROBO1$site == "01", ]
ROBO1 <- ROBO1[complete.cases(ROBO1), -3]
ROBO1$patient <- substr(ROBO1$sample, start = 1, stop = 12)

# plot distribution of expression levels
# mROBO1 <- melt(ROBO1, id = 'sample')
# p <- ggplot(mROBO1[mROBO1$variable == "ROBO1", ], aes(variable, value))
# p + geom_boxplot() + geom_jitter(width = 0.21)



# average sample by patient as group (one patient may have >=2 samples)
# ROBO1 <- aggregate(ROBO1[, 2], list(ROBO1$patient), mean)
# colnames(ROBO1) <- c("patient", "ROBO1")
#median ROBO1 level is 10.035
#25% 9.045; 75% 10.820
ROBO1$ROBO1_level <- "na"
ROBO1[ROBO1$ROBO1 >= 10.820, 'ROBO1_level'] <- "high"
ROBO1[ROBO1$ROBO1 < 10.820, 'ROBO1_level'] <- "low"
ROBO1 <- ROBO1[complete.cases(ROBO1), ]

#import clinical data
clinical <- read.delim("clinical.tsv", stringsAsFactors = FALSE)


#select clinical data
info <- clinical %>% select(submitter_id, days_to_death, days_to_last_follow_up, vital_status)
info$Died <- info$vital_status
info$Died <- gsub(info$Died, pattern = "dead", replacement = 1)
info$Died <- gsub(info$Died, pattern = "alive", replacement = 0)
info$Died <- as.numeric(info$Died)
info$survival.days <- ifelse(info$vital_status == "alive", info$days_to_last_follow_up, info$days_to_death)
info$survival.days <- as.numeric(info$survival.days)

#plot
plotdata <- merge(info, ROBO1, by.x = "submitter_id", by.y = "patient")
fit1 <- surv_fit(Surv(survival.days, Died) ~ ROBO1_level,data = plotdata)
ggsurvplot(fit1, data = plotdata, risk.table = TRUE)
surv_pvalue(fit1) 
