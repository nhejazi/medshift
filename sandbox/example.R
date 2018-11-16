library(data.table)
library(medshift)
library(mma)
data(weight_behavior)
missing <- unlist(apply(apply(weight_behavior,2, is.na), 2, which))
names(missing) <- NULL
missing <- unique(missing)

weight_data <- data.table(weight_behavior[-missing, ])
Y <- as.numeric(unlist(weight_data[, "overweigh"]))
A <- as.numeric(unlist(weight_data[, "snack"]))

