# ---------------------------------------------------------------------------- #
# Script 4B
# Author: Rani Powers
#
# Uses all features from pos & neg training examples to classify positive and 
# negative test examples. For a few training set sizes, store AUC values for all
# models to plot AUC and variability for each model.
# ---------------------------------------------------------------------------- #

cat('SCRIPT 4B: Train +/- (varying n), test +/- genes, seq features\n')

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

# Train on positive and negative genes, test on pos and neg genes
# For each training set size, store the AUC from each model

all.preds = list()
all.labels = list()

for (train.size in c(200, 100, 50, 10)){
  cat('Modeling with training size =', train.size, '\n')
  
  for (i in 1:10){
    cat('Iteration', i, '\n')
    
    # Divide data into train and test
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
    lcmix.perf$iteration = i
    
    if (exists('all.perf')){
      all.perf = rbind(all.perf, lcmix.perf)
    } else{
      all.perf = lcmix.perf
    }
    
    lcmix.preds = lcmix.res$preds
    
    # Add predictions to final list all.preds
    all.preds$new = lcmix.preds$semisup
    names(all.preds)[length(all.preds)] = paste0('Semisupervised ', train.size,
                                                 ' ', i)
    all.preds$new = lcmix.preds$unsup
    names(all.preds)[length(all.preds)] = paste0('Unsupervised ', train.size,
                                                 ' ', i)
    
    # Add labels to final list all.labels
    all.labels$new = labels[!is.train]
    names(all.labels)[length(all.labels)] = paste0('Semisupervised ', train.size,
                                                   ' ', i)
    all.labels$new = labels[!is.train]
    names(all.labels)[length(all.labels)] = paste0('Unsupervised ', train.size,
                                                   ' ', i)
    
    # Get J48 and LMT predictions
    rweka.train.data = cbind(divided.data$train,
                             sgd_ess = factor(labels[is.train]))
    J48.res = J48(sgd_ess~., data = rweka.train.data)
    LMT.res = LMT(sgd_ess~., data = rweka.train.data)
    J48.preds = predict(J48.res, divided.data$test, type = 'prob')
    LMT.preds = predict(LMT.res, divided.data$test, type = 'prob')
    
    # Add predictions to final list all.preds
    all.preds$new = J48.preds
    names(all.preds)[length(all.preds)] = paste0('C4.5 decision trees ', train.size,
                                                 ' ', i)
    all.preds$new = LMT.preds
    names(all.preds)[length(all.preds)] = paste0('Logistic model trees ', train.size,
                                                 ' ', i)
    
    # Add labels to final list all.labels
    all.labels$new = labels[!is.train]
    names(all.labels)[length(all.labels)] = paste0('C4.5 decision trees ', train.size,
                                                   ' ', i)
    all.labels$new = labels[!is.train]
    names(all.labels)[length(all.labels)] = paste0('Logistic model trees ', train.size,
                                                   ' ', i)
    
    # Get lasso predictions
    lasso.res = glmnet(as.matrix(divided.data$train), factor(labels[is.train]),
                       family = 'binomial')
    lasso.preds = predict(lasso.res, as.matrix(divided.data$test), 
                          type = 'response')[,length(lasso.res$lambda)]
    
    # Add predictions to final list all.preds
    all.preds$new = lasso.preds
    names(all.preds)[length(all.preds)] = paste0('Lasso ', train.size,
                                                 ' ', i)
    
    # Add labels to final list all.labels
    all.labels$new = labels[!is.train]
    names(all.labels)[length(all.labels)] = paste0('Lasso ', train.size,
                                                   ' ', i)
  }  
}

# Get AUC for each model
mroc = multiroc(all.preds, all.labels, quasi=F)
aucs = as.data.frame(do.call(rbind, lapply(mroc, function(l) l$auc)))
aucs$model = as.character(lapply(row.names(aucs), 
                                 function(x) strsplit(x, ' ')[[1]][1]))
aucs$train_size = as.numeric(lapply(row.names(aucs), function(x){ 
  s = strsplit(x, ' ')[[1]]
  return(s[length(s)-1])}))
aucs$iteration = as.numeric(lapply(row.names(aucs), function(x){ 
  s = strsplit(x, ' ')[[1]]
  return(s[length(s)])}))
names(aucs)[1] = 'AUC'

# Save AUC info
write.table(aucs[,c(2:4,1)], 'results/tables/4B-allFeatures_AUCs.txt', sep = '\t',
            row.names = F, quote = F)

# Save performance info
write.table(all.perf, 'results/tables/4B-allFeatures_pos_neg_lcmix_perf.txt',
            sep = '\t', row.names = F, quote = F)

# Plot AUCs
cols = brewer.pal(11, 'Spectral')
png('results/plots/4B-barplot_allFeatures_AUCs.png', height = 8, width = 15, 
    res = 300, units = 'in')
par(mfrow = c(1,4))
boxplot(AUC~model,data = aucs, subset = train_size == 200, 
        main = 'AUCs (train n = 200 +/- examples)', ylab = 'AUC', 
        ylim = c(.3, .7), col = cols[c(1,3,5,10,11)])
boxplot(AUC~model,data = aucs, subset = train_size == 100, 
        main = 'AUCs (train n = 100 +/- examples)', ylab = 'AUC', 
        ylim = c(.3, .7), col = cols[c(1,3,5,10,11)])
boxplot(AUC~model,data = aucs, subset = train_size == 50, 
        main = 'AUCs (train n = 50 +/- examples)', ylab = 'AUC', 
        ylim = c(.3, .7), col = cols[c(1,3,5,10,11)])
boxplot(AUC~model,data = aucs, subset = train_size == 10, 
        main = 'AUCs (train n = 10 +/- examples)', ylab = 'AUC', 
        ylim = c(.3, .7), col = cols[c(1,3,5,10,11)])
dev.off()