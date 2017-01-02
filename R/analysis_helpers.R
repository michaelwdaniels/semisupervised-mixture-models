# --------------------------------- WRAPPERS --------------------------------- #

# Functions to wrap the modeling and plotting portions of the analyses, which
# are repeated under lots of conditions (i.e. training set size)

divideYeastData <- function(feature.table, date.table, labels, 
                            train.date = NULL, 
                            train.class = c('positive', 'negative', 'both'),
                            num.train.examples = NULL){
  
  # Divide the yeast data set into a training and testing set based on desired
  # example class, date a positive example was discovered, or number of training
  # examples desired.
  
  # Format arguments
  train.class = train.class[1]
  if (!is.null(train.date) & train.class != 'positive'){ 
    stop('Dates only available when using positive training examples only')}
  if (train.class == 'both' & is.null(num.train.examples)){
    stop('If using both positive and negative training examples, must specify 
         the number of training examples to use')
  }
  
  # If using both pos and neg training examples, just split by the number of
  # training examples desired and we are done
  if (train.class == 'both'){
    is.train = row.names(feature.table) %in% sample(row.names(feature.table), 
                                                    num.train.examples)
    train.data = feature.table[is.train,]
    test.data = feature.table[!is.train,]
  } 
  
  # If using only pos training examples, first separate by label class
  if (train.class == 'positive'){
    train.data = feature.table[labels == 1,]
    test.data = feature.table[labels == 0,]
    
    # Then separate by date
    if (!is.null(train.date)){
      good.date = unique(row.names(date.table[date.table$year < train.date,]))
      is.good = row.names(train.data) %in% good.date
      test.data = rbind(test.data, train.data[!is.good,])
      train.data = train.data[is.good,]
    }
    
    # Then separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% sample(row.names(train.data), 
                                                        num.train.examples)
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }
  
  # If using only neg training examples, first separate by label class
  if (train.class == 'negative'){
    train.data = feature.table[labels == 0,]
    test.data = feature.table[labels == 1,]
    
    # The separate by number of training examples
    if (!is.null(num.train.examples)){
      if (num.train.examples < nrow(train.data)){
        is.train = row.names(train.data) %in% sample(row.names(train.data), 
                                                     num.train.examples)
        test.data = rbind(test.data, train.data[!is.train,])
        train.data = train.data[is.train,]
      } else{
        stop('Number of training examples desired is larger than training set')
      }
    }
    is.train = row.names(feature.table) %in% row.names(train.data)
  }

  # Return divided data
  results = list(feature.table = feature.table,
                 is.train = is.train,
                 train = train.data,
                 test = test.data)
  return(results)
}

getLCMIXmodels <- function(feature.table, is.train, test.labels, K0 = 2, TW = 10){
  
  BINARY_FEATURES2 = BINARY_FEATURES[BINARY_FEATURES %in% names(feature.table)]
  INTEGER_FEATURES2 = INTEGER_FEATURES[INTEGER_FEATURES %in% names(feature.table)]
  OPEN_FEATURES2 = OPEN_FEATURES[OPEN_FEATURES %in% names(feature.table)]
  CLOSED_FEATURES2 = CLOSED_FEATURES[CLOSED_FEATURES %in% names(feature.table)]
  
  # Make feature tables
  X = list(
    bern = do.call(cbind, feature.table[,BINARY_FEATURES2]),
    nbin = do.call(cbind, feature.table[,INTEGER_FEATURES2]),
    norm = cbind(
      do.call(cbind, feature.table[,OPEN_FEATURES2]),
      log(do.call(cbind, feature.table[,CLOSED_FEATURES2]))
    )
  )
  
  # Marginal model selection
  message(notification("marginal model selection", 1))
  margsel = msapply(namedList(2,3,4,5), function(K)
    mapply(X, names(X), FUN=function(x, fam)
      iclbic(mixmod(x, K=K, family=fam), map=T)))
  K = apply(margsel, 1, which.max) + 1
  
  # Joint model fitting
  message(notification("fitting joint models", 1))
  train = as.numeric(is.train)
  fm1 = mcparallel(mdmixmod(X, K=K, K0=K0, family=names(X),
                            train=train, train.weight=TW)) 
  fm2 = mcparallel(mdmixmod(X, K=K, K0=K0, family=names(X)))
  fits = mccollect(list(fm1, fm2))
  names(fits) = c("semisup", "unsup")
  
  # Make predictions for semisup and unsup models
  preds = lapply(fits, function(x) predict(x, type="prob")[!is.train, 1])
  
  # Report performance
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