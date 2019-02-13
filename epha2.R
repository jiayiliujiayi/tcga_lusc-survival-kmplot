setwd("/Users/JiayiLiu/Downloads/ngsanalysis/")

library(dplyr);library(stringr)
#import EPHA2 expression level, add patient barcode
originalEPHA2 <- read.delim("EPHA2.tsv", stringsAsFactors=FALSE)

EPHA2 <- originalEPHA2
EPHA2 <- EPHA2[complete.cases(EPHA2), ]
EPHA2$patient <- substr(EPHA2$sample, start = 1, stop = 12)
#average sample by patient as group (one patient may have >=2 samples)
EPHA2 <- aggregate(EPHA2[, 2], list(EPHA2$patient), mean)
colnames(EPHA2) <- c("patient", "EPHA2")
#median EPHA2 level is 10.6050
#25% 9.9165; 75% 11.2725
EPHA2$EPHA2_level <- "na"
EPHA2[EPHA2$EPHA2 >= 11.2725, 3] <- "high"
EPHA2[EPHA2$EPHA2 < 11.2725, 3] <- "low"
EPHA2 <- EPHA2[complete.cases(EPHA2), ]

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
plotdata <- merge(info, EPHA2, by.x = "submitter_id", by.y = "patient")
fit1 <- surv_fit(Surv(survival.days, Died) ~ EPHA2_level,data = plotdata)
ggsurvplot(fit1, data = plotdata, risk.table = TRUE)
surv_pvalue(fit1) 
