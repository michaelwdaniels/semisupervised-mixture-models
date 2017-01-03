# ---------------------------------------------------------------------------- #
# Script 3B
# Author: Rani Powers
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

# Train on positive and negative genes, test on pos and neg genes

all.preds = list()
all.labels = list()
for (train.size in c(200, 100, 50, 10)){
  
  # Divide data into train and test
  cat('Modeling with training size =', train.size, '\n')
  divided.data = divideYeastData(feature.table, date.table, labels,
                                 train.class = 'both', 
                                 num.train.examples = train.size)
  is.train = divided.data$is.train
  
  # Make sure both classes are represented in training set
  cat(table(labels[is.train]), '\n')
  while (length(table(labels[is.train])) != 2 | 
         any(as.numeric(table(labels[is.train])) < 2)){
    divided.data = divideYeastData(feature.table, date.table, labels,
                                   train.class = 'both', 
                                   num.train.examples = train.size)
    is.train = divided.data$is.train
    cat(table(labels[is.train]), '\n')
  }
  
  # Get lcmix performance and predictions
  lcmix.res = getLCMIXmodels(feature.table, is.train, labels[!is.train])
  lcmix.perf = lcmix.res$performance
  lcmix.perf$train_size = train.size
  
  if (exists('all.perf')){
    all.perf = rbind(all.perf, lcmix.perf)
  } else{
    all.perf = lcmix.perf
  }
  
  lcmix.preds = lcmix.res$preds
  
  # Add predictions to final list all.preds
  all.preds$new = lcmix.preds$semisup
  names(all.preds)[length(all.preds)] = paste0('Semisupervised (train n = ', 
                                               train.size, ')')
  all.preds$new = lcmix.preds$unsup
  names(all.preds)[length(all.preds)] = paste0('Unsupervised (train n = ', 
                                               train.size, ')')
  
  # Add labels to final list all.labels
  all.labels$new = labels[!is.train]
  names(all.labels)[length(all.labels)] = paste0('Semisupervised (train n = ', 
                                                 train.size, ')')
  all.labels$new = labels[!is.train]
  names(all.labels)[length(all.labels)] = paste0('Unsupervised (train n = ', 
                                                 train.size, ')')
  
  # Get J48 and LMT predictions
  rweka.train.data = cbind(divided.data$train,
                           sgd_ess = factor(labels[is.train]))
  J48.res = J48(sgd_ess~., data = rweka.train.data)
  LMT.res = LMT(sgd_ess~., data = rweka.train.data)
  J48.preds = predict(J48.res, divided.data$test, type = 'prob')
  LMT.preds = predict(LMT.res, divided.data$test, type = 'prob')
  
  # Add predictions to final list all.preds
  all.preds$new = J48.preds
  names(all.preds)[length(all.preds)] = paste0('C4.5 decision trees (train n = ', 
                                               train.size, ')')
  all.preds$new = LMT.preds
  names(all.preds)[length(all.preds)] = paste0('Logistic model trees (train n = ', 
                                               train.size, ')')
  
  # Add labels to final list all.labels
  all.labels$new = labels[!is.train]
  names(all.labels)[length(all.labels)] = paste0('C4.5 decision trees (train n = ', 
                                                 train.size, ')')
  all.labels$new = labels[!is.train]
  names(all.labels)[length(all.labels)] = paste0('Logistic model trees (train n = ', 
                                                 train.size, ')')
  
  # Get lasso predictions
  lasso.res = glmnet(as.matrix(divided.data$train), factor(labels[is.train]),
                     family = 'binomial')
  lasso.preds = predict(lasso.res, as.matrix(divided.data$test), 
                        type = 'response')[,length(lasso.res$lambda)]
  
  # Add predictions to final list all.preds
  all.preds$new = lasso.preds
  names(all.preds)[length(all.preds)] = paste0('Lasso (train n = ', 
                                               train.size, ')')
  
  # Add labels to final list all.labels
  all.labels$new = labels[!is.train]
  names(all.labels)[length(all.labels)] = paste0('Lasso (train n = ', 
                                                 train.size, ')')
  
}

# Save performance info
write.table(all.perf, 'results/tables/3B-allFeatures_pos_neg_lcmix_perf.txt',
            sep = '\t', row.names = F, quote = F)

# Plot all lines together
mroc = multiroc(all.preds, all.labels, quasi=F)
cols = brewer.pal(11, 'Spectral')
png('results/plots/3B-ROC_allFeatures_pos_neg.png',
    height = 8, width = 8, res = 300, units = 'in')
plot(mroc, grid=TRUE, col=rep(c('red', 'black', 'blue', 'green', 'purple'), 4),
     cut=c(0.5, 0.8, 0.99), 
     main = 'Train: essential and non-essential genes (varying n)\nTest: essential and non-essential genes\nData: 22 Seringhaus features',
     cex.legend = .5)
dev.off()

# Plot only lcmix lines
mroc.lcmix = multiroc(all.preds[grepl('supervised', names(all.preds))], 
                      all.labels[grepl('supervised', names(all.labels))], quasi=F)
png('results/plots/3B-ROC_allFeatures_pos_neg_lcmix.png',
    height = 8, width = 8, res = 300, units = 'in')
plot(mroc.lcmix, grid=TRUE, col=rep(c('red', 'black'), 4),
     cut=c(0.5, 0.8, 0.99), 
     main = 'Train: essential and non-essential genes (varying n)\nTest: essential and non-essential genes\nData: 22 Seringhaus features')
dev.off()
