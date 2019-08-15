# This script creates supervised models of response using the microbial data in the 
# melanoma data sets (Gopalakrishnan et al., Matson et al., and Frankel et al.) as
# training data. The models are then use to predict response in RCC and NSCLC data
# (Routy et al.).

# ==== LOAD DEPENDENCIES ====

library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("Rtsne")
library("dplyr")
library("data.table")
library("tibble")
library("randomForest")
library("EnvStats")
library("reshape2")
library("pROC")
library("xgboost")
library("svglite")
library("fastDummies")
## if not installed, quickly add it as follows:
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("RColorBrewer", "pheatmap", "ggplot2", "Rtsne"))

# ==== IMPORT DATA ====

# Abundance table for each OTU
OTU = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/OTU.rds")
# Clinical metadata
clin = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")
# Pyholgenetic metadata
phylo_meta = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/phylo_meta.rds")

# ==== DATA PRE-PROCESSING ====

# Convert to relative abundance
OTU = apply(OTU, 2, function(x){x/sum(x)})
# Filter out microbes that are not identified across at least 5 samples
OTU_filtered = OTU[apply(OTU, 1, function(x) sum(x > 0)) > 5,]
# log transform the data
OTU_filtered = log10(OTU_filtered + 1e-8)
# Sort according to most abundant OTUs
OTU_filtered = OTU_filtered[order(rowSums(OTU_filtered),decreasing = TRUE),]

# Split into training and test data
clin_train = clin %>% filter(Study != "Routy_et_al")
OTU_train = as.data.frame(t(OTU_filtered[,colnames(OTU_filtered) %in% clin_train$Patient_id]))
clin_test = clin %>% filter(Study == "Routy_et_al")
OTU_test = as.data.frame(t(OTU_filtered[,colnames(OTU_filtered) %in% clin_test$Patient_id]))

#Concatenate Treatment, Study and Cancer_type information into data frame
OTU_train = cbind(OTU_train,
                  Treatment = clin_train$Treatment,
                  Study = clin_train$Study,
                  Cancer_type = clin_train$Cancer_type)
OTU_test = cbind(OTU_test,
                  Treatment = clin_test$Treatment,
                  Study = clin_test$Study,
                  Cancer_type = clin_test$Cancer_type)

# ==== RANDOM FOREST ====
set.seed(8)
# Construct random forest model for training data
rf.OTU <- randomForest(x=OTU_train,y=clin_train$Response, ntree = 10000 ,importance = TRUE)
# View confusion matrix
print(rf.OTU$confusion)
# View the most important (Mean decrease of GINI impurity) features for prediction
rf.OTU.importance <- data.frame(importance(rf.OTU, type = 2)) %>%
  rownames_to_column('mOTU') %>%
  arrange(desc(MeanDecreaseGini))
print(head(rf.OTU.importance))

# Let's plot the top predictiors
nPlot = 12
# Select top predictors from OTU matrix and format for ggplot
plotData = as.matrix(OTU_train[,rf.OTU.importance$mOTU[1:nPlot]])
plotData = melt(plotData)
plotData = cbind(plotData, Response = gsub(".*_R","R",gsub(".*_NR","NR",plotData$Var1)))
colnames(plotData) = c("Patient","mOTU","Abundance","Response")
p_RFpredictors = ggplot(data = plotData,aes(x=Response,y=Abundance,color=Response)) +
  facet_wrap(facets = "mOTU",ncol = 4, scales = "free") +
  geom_violin() +
  geom_jitter()
p_RFpredictors

# Generate ROC-curve for training data
rf.train.roc = roc(as.factor(clin_train$Response),rf.OTU$votes[,2])
plot(rf.train.roc)
auc(rf.train.roc)

# That's pretty bad, lets try vaying number of trees
n_tree = c(50, 100, 500, 1000, 2000, 10000, 20000)
plotData = data.frame(n_tree, AUC = 1:length(n_tree))
for (i in 1:length(n_tree)) {
  set.seed(8)
  print(paste("number of trees: ",n_tree[i]))
  rf.OTU <- randomForest(x=OTU_train,y=as.factor(clin_train$Response), ntree = n_tree[i] ,importance = TRUE)
  rf.train.roc = roc(as.factor(clin_train$Response),rf.OTU$votes[,2])
  plotData$AUC[i] = auc(rf.train.roc)
}

p_RFperf = ggplot(data = plotData, aes(x = n_tree, y = AUC)) +
  geom_line() +
  geom_point() +
  xlab("no. of trees") +
  ylab("AUC") +
  ggtitle("Random forest on OTU matrix") +
  theme(plot.title = element_text(face="bold")) +
  theme(text = element_text(family = "Helvetica"))
p_RFperf

# Save plots
ggsave(filename = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Figures/RFpredictors.svg",
       plot = p_RFpredictors,
       width = 9, 
       height = 8, 
       pointsize = 10)
ggsave(filename = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Figures/RFperformance.svg",
       plot = p_RFperf,
       width = 4, 
       height = 4, 
       pointsize = 10)

# ==== GRADIENT BOOSTING TREES ====

# Let's reformat our data for XGboost
# We will use one-hot encoding for factor variables
OTU_train_xgb = dummy_cols(.data = OTU_train, select_columns = names(Filter(is.factor, OTU_train))) %>%
  select(-one_of(names(Filter(is.factor, OTU_train)))) # Drop all factor columns

params_xgb = list(eta = 0.3,
                  max_depth =6)
set.seed(8)
xgb.OTU = xgb.cv(data = as.matrix(OTU_train_xgb),
                 label = as.numeric(clin_train$Response)-1,
                 nthread = 4,
                 nrounds = 3,
                 nfold = 5,
                 num_parallel_tree = 1000,
                 objective = "binary:logistic",
                 metrics = list("rmse","auc"))
print(xgb.OTU)

xgb.plot.tree(model = xgb.OTU)
plot(rf.probTrain.performance)

CairoWin()
barplot(rf.iris.importance, beside = TRUE,
        names.arg = rownames(rf.iris.importance),
        col = c("red","cornflowerblue","green2","mediumvioletred"),
        ylab = "Mean decrease of GINI impurity")