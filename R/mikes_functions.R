
getFeatureTable <- function(feature=c('seq','all')){
  load('data/Rdata/cerevisiae_preprocessed.RData')
  
  # Format feature.table, date.table, and labels for divideYeastData()
  
  feature.table = as.data.frame(preproc$features)
  
  if (feature == 'seq'){
    feature.table = feature.table[,SEQUENCE_DERIVED_FEATURES]
  } 
  if (feature == 'all'){
    feature.table = feature.table
  } 
  row.names(feature.table) = preproc$ensembl
  return(feature.table)
}

GetTrainTestSet <- function(feature.table=feature.table,labels=labels,p=NULL,contam=NULL,design=c('balanced','unbalanced'),seed=seed){
  # p = percent training level, contam = percent contamination, design = balanced or unbalanced testing, seed = seed
  # First identify row numbers for positive and negative labels
  labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
  labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
  
  # Second determine which positive labels will be used in the training set and how many
  number.pos <- floor(p*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
  idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
  inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
  
  # Third determine number of negative labels and the subset number of positive labels to be the contamination
  number.neg <- floor(p*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
  number.contam <- floor(contam*number.neg) # among the number of negatives how many will be contamination based on contamination percentages
  contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
  contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
  
  # Fourth determine true negative labels and assemble full training set
  idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
  inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
  inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
  inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
  training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
  
  # Training Set must not contain a column of zeros or featured variable will be treated as a constant
  repeat{
    if(product(apply(training[,1:7],2,sum))!=0) {
      print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
      break
    } else{ 
      seed <- seed+j*i
      set.seed(seed)
      idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
      inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
      contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
      contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
      idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
      inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
      inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
      training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
    }
  }   
  
  # Construct test set by sampling equal numbers from remaining labels
  # The test set should contain the same number of positive labels as is in the training set
  # and the same number of negative labels as the combined number of true negatives and contaminates in the training set
 if(design=='balanced'){ 
  # Fifth collect all row numbers for remaining positives and negatives
  labels.pos.test <- labels.pos[-c(idx.pos,contam.idx)] # collection of row numbers from labels which are positive but not used in training set (includes true positives and contamination)
  labels.neg.test <- labels.neg[-idx.neg] # collection of row numbers from labels which are negative
  
  # Sixth sample from test labels with corresponding training sizes (balanced approach)
  idx.pos.test <- sample(1:length(labels.pos.test), size=number.pos[i], replace=FALSE) # sample row numbers from remaining positive labels
  idx.neg.test <- sample(1:length(labels.neg.test), size=number.neg[i], replace=FALSE) # sample row numbers from remaining negative labels
  idx.test <- c(labels.pos.test[idx.pos.test],labels.neg.test[idx.neg.test]) # assemble row numbers of full testing set
  testing <- feature.table[idx.test,] # Testing set constructed using row number ids within feature table
 }
  else{testing <- feature.table[,-inTrain]}
}

getPerf <- function(preds=NULL,test.labels=NULL){
reports = c("precis", "recall", "f1meas")
suppressWarnings(rept1 <- as.data.frame(sapply(preds, perfmeas, 
                                               labels=test.labels)[reports,]))
rept1$report = reports
rept1$FDR = 0
suppressWarnings(rept2 <- as.data.frame(sapply(preds, perfmeas, lfdr=.20, 
                                               labels=test.labels)[reports,]))
rept2$report = reports
rept2$FDR = .2
suppressWarnings(rept3 <- as.data.frame(sapply(preds, perfmeas, lfdr=.01, 
                                               labels=test.labels)[reports,]))
rept3$report = reports
rept3$FDR = .01
combined.report = rbind(melt(rept1, id.vars = c('report', 'FDR'), 
                             measure.vars = c('semisup', 'unsup')),
                        melt(rept2, id.vars = c('report', 'FDR'), 
                             measure.vars = c('semisup', 'unsup')),
                        melt(rept3, id.vars = c('report', 'FDR'), 
                             measure.vars = c('semisup', 'unsup')))

# Return prediction probabilities and performance report
return(list(preds = preds, performance = combined.report))
}

