library(lubridate)
library(tidyverse)

metadata <- readxl::read_xlsx("metadata.xlsx")
metadata$Age <- interval(dmy(metadata$DOB), dmy(metadata$StoolCollectionDate)) / years(1)
metadata$AgeGroup <- ifelse(metadata$Age <=3, "Young", "Old")
metadata$Race <- replace(metadata$Race, metadata$Race %in% c("Eurasian", "Caucasian"), "Others")
metadata$FeedingMethod <- replace(metadata$FeedingMethod,is.na(metadata$FeedingMethod), "None")
write_tsv(metadata, file='metadata.tsv')
