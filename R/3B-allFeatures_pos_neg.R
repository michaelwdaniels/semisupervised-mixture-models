# ---------------------------------------------------------------------------- #
# Script 3B
# Authors: Rani Powers and Daniel Dvorkin
#
# Uses all features from pos & neg training
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 3B: Train +/- (varying n), test +/- genes, all features\n')

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')

# Load preproc - preprocessed data from Script 0
load('data/Rdata/cerevisiae_preprocessed.RData')

# Format feature.table, date.table, and labels for divideYeastData()
feature.table = as.data.frame(preproc$features)
row.names(feature.table) = preproc$ensembl

date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #
# --------------------------- GET LCMIX PREDICTIONS -------------------------- #

# Train on positive genes, test on pos and neg genes

all.preds = list()
all.labels = list()
for (train.size in c(200, 100, 50, 30, 20, 10, 5)){
  cat('Modeling with training size =', train.size, '\n')
  divided.data = divideYeastData(feature.table, date.table, labels,
                                 train.class = 'both', 
                                 num.train.examples = train.size)
  is.train = divided.data$is.train
  lcmix.preds = getLCMIXmodels(feature.table, is.train, labels[!is.train])$preds
  all.preds$new = lcmix.preds$semisup
  names(all.preds)[length(all.preds)] = paste0('Semisupervised (train n = ', 
                                               train.size, ')')
  all.preds$new = lcmix.preds$unsup
  names(all.preds)[length(all.preds)] = paste0('Unsupervised (train n = ', 
                                               train.size, ')')
  all.labels$new = labels[!divided.data$is.train]
  names(all.labels)[length(all.labels)] = paste0('Semisupervised (train n = ', 
                                                 train.size, ')')
  all.labels$new = labels[!divided.data$is.train]
  names(all.labels)[length(all.labels)] = paste0('Unsupervised (train n = ', 
                                                 train.size, ')')
}


# Plot all lines together
mroc = multiroc(all.preds, all.labels, quasi=F)
cols = brewer.pal(11, 'Spectral')
png('results/plots/ROC_allFeatures_pos_neg.png',
    height = 8, width = 8, res = 300, units = 'in')
plot(mroc, grid=TRUE, col=rep(c('red', 'black'), 4),
     cut=c(0.5, 0.8, 0.99), 
     main = 'Train: essential and non-essential genes (varying n)\nTest: essential and non-essential genes\nData: 22 Seringhaus features')
dev.off()

rm(list = ls())