install.packages("powerMediation")
library(powerMediation)

minEffect.VSMc.cox()

qnorm(0.95)

setwd("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/")

val_clin = readRDS("./Metadata/Processed_metadata/val_clin.rds")




delta=0
# The idea is that statistically significant differences
# between the proportions may not be of interest unless
# the difference is greater than a threshold, δ. This is
# particularly popular in clinical studies, where the 
# margin is chosen based on clinical judgement and subject-
# domain knowledge.
pA = sum(val_clin$Response == "R")/nrow(val_clin)
alpha = 0.05
beta = 0.20
n = 27
pE = sum(!is.na(val_clin$months_to_progression))/n


logHazard = delta - ((qnorm(1-alpha)+qnorm(1-beta/2))/sqrt(n*pA*(1-pA)*pE))
print(paste0("log hazard-ratio: ", logHazard))
print(paste0("hazard-ratio: ", exp(logHazard)))


