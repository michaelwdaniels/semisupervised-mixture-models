# ---------------------------------------------------------------------------- #
# Script 1A
# Author: Rani Powers
#
# Uses the 14 Seringhaus sequence-derived features from the pre-2002 essential 
# genes to predict which post-2002 genes are essential.
# ---------------------------------------------------------------------------- #

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')

cat('SCRIPT 1A: Train pre-2002 (+ only), test 2002+ (+/- genes), sequence-derived features\n')

# Load preproc - preprocessed data from Script 0
load('data/Rdata/cerevisiae_preprocessed.RData')

# Format feature.table, date.table, and labels for divideYeastData()
feature.table = as.data.frame(preproc$features)
feature.table = feature.table[,SEQUENCE_DERIVED_FEATURES]
row.names(feature.table) = preproc$ensembl

date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #

# Train on the pre-2002 pos genes, test on the 2002+ pos and neg genes
divided.data = divideYeastData(feature.table, date.table, labels,
                               train.date = 2002, train.class = 'positive')


# --------------------------- GET LCMIX PREDICTIONS -------------------------- #

cat('Making predictions with lcmix\n')
is.train = divided.data$is.train
lcmix.preds = getLCMIXmodels(feature.table, is.train, labels[!is.train])

mroc = multiroc(lcmix.preds$preds, labels[!divided.data$is.train], quasi=F)
names(mroc) = c('Semisupervised' ,'Unsupervised')
png('results/plots/1A-ROC_seqFeatures_prepost2002.png',
    height = 8, width = 8, res = 300, units = 'in')
plot(mroc, grid=TRUE, col=c("red", "black"),
     cut=c(0.5, 0.8, 0.99), 
     main = 'Train: pre-2002 essential genes (n = 64)\nTest: 2002+ essential and non-essential genes (n = 3,436)\nData: 14 sequence-derived features')
dev.off()

# --------------------- RWEKA MODELING ESSENTIAL GENES ----------------------- #

#cat('Making predictions with J48 (C4.5 decision trees)\n')

#rweka.train.data = cbind(divided.data$train, 
#                         sgd_ess = factor(labels[divided.data$is.train]))

# J48, LMT and M5P cannot handle unary class labels
#j48.res = J48(sgd_ess~., data = rweka.train.data)
#LMT.res = LMT(sgd_ess~., data = rweka.train.data)
#M5P.res = M5P(sgd_ess~., data = rweka.train.data)

#library(e1071)
#svm.model = svm(divided.data$train,
#               y=NULL,
#               type='one-classification',
#               nu=0.5,
#               scale=F,
#               kernel="radial")
#svm.pred = predict(svm.model,divided.data$test)

