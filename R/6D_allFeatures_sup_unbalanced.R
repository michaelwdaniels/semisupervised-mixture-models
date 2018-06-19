setwd("~/Desktop/semisupervised-mixture-models-master")
load("saved_images/6D_boxplots_unbalanced")
load('data/Rdata/cerevisiae_preprocessed.RData')
# ---------------------------------------------------------------------------- #
# Script 6D all unbalanced
# Author: Mike Daniels
#
# Uses all features from only positive training for semi-supervised method
# compares to supervised methods (SVM, Decision Tree, LASSO)
# Negatives in training set will be contaminated with various percentages of positives to reflect a real world application
# ---------------------------------------------------------------------------- #

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')
source('R/mikes_functions.R')
library(randomForest)
library(caret)
library(penalized)
library(e1071)
library(glmnet)
library(ROCR)
library(gridExtra)
library(grid)
library(MAd)
library(sm)
library(RSSL)

feature.table <- getFeatureTable(feature='all')
dim(feature.table)
names(feature.table)
str(feature.table)
head(feature.table)
apply(feature.table[,c(2,7,9,10,13,17,20,21,22)],2,sum)/3500
#########################
# remove integer variables with less than 5% information
idx <- which(apply(feature.table[,c(2,7,9,10,13,17,20,21,22)],2,sum)/3500<0.05)
# removed from analysis due to low content -  vacuole, in_how_many_of_5_proks_blast, intron 
idx <- c(2,7,9,10,13,17,20,21,22)[idx]
names(feature.table)[idx]
feature.table <- feature.table[,-idx]
########################
names(feature.table) # 3500 19
date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# Perform iterative calculations using supervised and semisupervised
n1 <- 100  # Number of iterations
p <- c(seq(0.01,0.1,0.005)) # percentage training set values
l <- length(p)
n <- p*length(labels)
contam.perc <- seq(0,0.9,0.1)

AUC.svm <- matrix(0,ncol=n1,nrow=l) # collects AUC values for SVM analysis
AUC.J48 <- matrix(0,ncol=n1,nrow=l) # collects AUC values for decision tree analysis
AUC.sup <- matrix(0,ncol=n1,nrow=l) # collects AUC values for LASSO analysis
AUC.semi <- matrix(0,ncol=n1,nrow=l) # collects AUC values for semi-supervised analysis
label.svm <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for SVM for contamination confirmation 
label.tree <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for decision tree for contamination confirmation
label.lasso <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for LASSO for contamination confirmation
# Contamination does not affect semi-supervised analysis because only positive labels are used in the training set 
rownames(AUC.J48) <- rownames(AUC.svm) <- rownames(AUC.sup) <- rownames(AUC.semi) <- p # row names are training percent levels, columns will be iterations

constant=sum(labels) #769 total number of positive labels in yeast data

sup.prec.med <- sup.prec.pp50 <- sup.prec.pp80 <- sup.prec.pp95 <- sup.prec.midrange <-
  sup.rec.med <- sup.rec.pp50 <- sup.rec.pp80 <- sup.rec.pp95 <- sup.rec.midrange <- 
  sup.f.med <- sup.f.pp50 <- sup.f.pp80 <- sup.f.pp95 <- sup.f.midrange <- 
  sup.pos.med <- sup.pos.pp50 <- sup.pos.pp80 <- sup.pos.pp95 <- sup.pos.midrange <- 
  semi.prec.med <- semi.prec.pp50 <- semi.prec.pp80 <- semi.prec.pp95 <- semi.prec.midrange <-
  semi.rec.med <- semi.rec.pp50 <- semi.rec.pp80 <- semi.rec.pp95 <- semi.rec.midrange <- 
  semi.f.med <- semi.f.pp50 <- semi.f.pp80 <- semi.f.pp95 <- semi.f.midrange <- 
  semi.pos.med <- semi.pos.pp50 <- semi.pos.pp80 <- semi.pos.pp95 <- semi.pos.midrange <- matrix(0,nrow = l,ncol = n1) 

AUC.list <-  sup.perf <- sup.mean <- sup.median <- sup.var <- sup.mad <- forest.mean <- 
  forest.median <- forest.var <- forest.mad <- svm.mean <- svm.median <- svm.var <- svm.mad <- all.pred.svm <- 
  all.pred.semi <- all.pred.lasso <- all.pred.forest <- list()

for(g in 1:length(contam.perc)){#length(contam.perc)

cat(paste0('Iteration g: ',g,' of ',length(contam.perc),'\n'))

    for(j in 1:n1){#n1

    cat(paste0('Iteration g: ',g,' of ',length(contam.perc),', j: ',j,' of ',n1,'\n'))
    
      # Create lists to be used in analysis
    all.preds.svm <- list() # list of predictions
    all.labels.svm <- list() # list of labels being predicted
    all.preds.J48 <- list() # list of predictions
    all.labels.J48 <- list() # list of labels being predicted
    all.preds.sup.p <- list() # list of predictions
    all.labels.sup <- list() # list of labels being predicted
   if(g<2) {
    all.preds.semi <- list()
    all.labels.semi <- list()
    all.perfs.semi <- list()}
    n.pos = list() # list of number of positive labels in each iteration
    number.pos <- vector() # list of number of positive labels in training set
    number.neg <- vector() # list of total number of negative labels in training set
    number.contam <- vector() # list of number of contaminates in training set
    inTrain.neg.only.length <- vector() # list of number of true negative labels in training set
    
    # Within each iteration, perform analysis for each training percent level (l is the number of training sizes being analyzed)
    for(i in 1:l){#l
  
      cat(paste0('Iteration  g: ',g,' of ',length(contam.perc),', j: ',j,' of ',n1,', i: ',i,' of ',l,'\n'))
      
      seed <- 210 + g*i*j
      seed1 <- seed
      set.seed(seed)
      # First identify row numbers for positive and negative labels
      labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
      labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
      
      # Second determine which positive labels will be used in the training set and how many
      number.pos[i] <- floor(p[i]*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
      idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
      n.pos$new <- number.pos[i] # collect the number of positive labels in a list for contamination calculations
      names(n.pos)[length(n.pos)] = paste0('train n = ', floor(p[i]*length(labels))) # each assignment must be labeled in order to start a new list entry
      inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
      
      # Third determine number of negative labels and the subset number of positive labels to be the contamination
      number.neg[i] <- floor(p[i]*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
      number.contam[i] <- floor(contam.perc[g]*number.neg[i]) # among the number of negatives how many will be contamination based on contamination percentages
      contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam[i],replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
      contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
      labels.contam <- labels
      labels.contam[contam] <- 0
      
      # Fourth determine true negative labels and assemble full training set
      idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
      inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
      inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
      inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
      training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
      testing <- feature.table[-inTrain,] # Testing set constructed using row number ids within feature table
      
      repeat{
        if(product(apply(training[,c(2,7,9,12)],2,sum))!=0 & product(apply(testing[,c(2,7,9,12)],2,sum))!=0) {
          print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
          break
        } else{ 
          seed <- 210 + seed
          set.seed(seed)
          idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
          inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
          contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
          contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
          labels.contam <- labels
          labels.contam[contam] <- 0
          idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
          inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
          inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
          training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
          testing <- feature.table[-inTrain,] # Testing set constructed using row number ids within feature table
        }
      }  
      
      # SVM modeling 
     
      svm.model = svm(training, y=as.factor(labels.contam[inTrain]), kernel="radial",probability = TRUE) # radial kernel showed the best performance from initial trials
      all.preds.svm$new = predict(svm.model,testing,type="C-classification",probability=TRUE) # predictions of svm model on test set
      names(all.preds.svm)[length(all.preds.svm)] = paste0('SVM (train n = ', floor(p[i]*length(labels)), ')') # assign names to predictions in a list
      all.labels.svm$new = labels[-inTrain] # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.svm)[length(all.labels.svm)] = paste0('SVM (train n = ', floor(p[i]*length(labels)), ')') # assign names to labels in a list
      
      # Random Forest
      training.tree <- cbind(labels=as.factor(labels.contam[inTrain]),training)
      J48.model = randomForest(labels~., data = training.tree)
      predict.J48 <- predict(J48.model,testing,type = "class") # when type = "class", a factor vector is returned. When type = "prob", a matrix of confidence values is returned (one column per class).
      all.preds.J48.p <- predict(J48.model,testing,type="prob")[,1]
      labels.J48 <- labels[-inTrain]
      
      all.preds.J48$new = predict.J48 # predictions of tree model on test set
      names(all.preds.J48)[length(all.preds.J48)] = paste0('Tree (train n = ', floor(p[i]*length(labels)), ')') # assign names to predictions in a list
      all.labels.J48$new = labels.J48 # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.J48)[length(all.labels.J48)] = paste0('Tree (train n = ', floor(p[i]*length(labels)), ')') # assign names to labels in a list
      
      # LASSO modeling 
      cv.lasso <- cv.glmnet(as.matrix(training), as.double(labels.contam[inTrain]), family='binomial', alpha=1, standardize=TRUE, type.measure='auc')
      
      fit <- glmnet(as.matrix(training), as.double(labels.contam[inTrain]), family='binomial',alpha = 1,lambda=cv.lasso$lambda.min)    
      
      all.preds.sup.p$new = predict(fit,as.matrix(testing),type="response",probability=TRUE) # predictions of sup model on test set
      
      names(all.preds.sup.p)[length(all.preds.sup.p)] = paste0('LASSO (train n = ', floor(p[i]*length(labels)), ')') # assign names to predictions in a list
      all.labels.sup$new = labels[-inTrain] # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.sup)[length(all.labels.sup)] = paste0('LASSO (train n = ', floor(p[i]*length(labels)), ')') # assign names to labels in a list

      
     if(g<2){
        # Semi-sup modeling
        train.size <- floor(p[i]*length(labels)) ###### training sets sizes were assumming only 769 total genes
        
        cat("Modeling with training size =", train.size, "\n")
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = "positive",
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = "all")
        
        repeat{
          if(product(apply(divided.data$train[,c(2,7,9,13)],2,sum))!=0 &
             product(apply(divided.data$test[,c(2,7,9,13)],2,sum))!=0 ) {
            print(paste0("Seed1 = ",seed1," doesnt contain a column of zeros."))
            break
          } else{ 
            seed1 <- 210 + seed1 
            set.seed(seed1)
            divided.data = divideYeastData(feature.table, date.table, labels,
                                           train.class = "positive",
                                           num.train.examples = train.size, seed = seed1, constant = constant, negatives = "all")
            }
        }   
        
        is.train = divided.data$is.train
        
        lcmix.out = getLCMIXmodels.semi.only(feature.table, is.train, labels[!is.train])
        
        all.preds.semi$new = lcmix.out$preds$semisup
        names(all.preds.semi)[length(all.preds.semi)] = paste0("Semisupervised (train n = ", 
                                                               train.size, ")")
        all.labels.semi$new = labels[!divided.data$is.train]
        names(all.labels.semi)[length(all.labels.semi)] = paste0("Semisupervised (train n = ", 
                                                                 train.size, ")")
         }
      
   
   if(g==1|g==3|g==6){
     
     sup.pred1 <- prediction(all.preds.sup.p[[i]],all.labels.sup[[i]]) # predicted response vectors
     # sup.pred <- slot(sup.pred1,'predictions')[[1]]/max(slot(sup.pred1,'predictions')[[1]])
     if(max(slot(sup.pred1,'predictions')[[1]])-min(slot(sup.pred1,'predictions')[[1]])>0){
     sup.pred <- (slot(sup.pred1,'predictions')[[1]]-min(slot(sup.pred1,'predictions')[[1]]))/(max(slot(sup.pred1,'predictions')[[1]])-min(slot(sup.pred1,'predictions')[[1]]))
     
     idx.sup <- order(sup.pred)
     sup.pred.sorted <- sup.pred[idx.sup]
     sup.lab.sorted <- as.numeric(slot(sup.pred1,'labels')[[1]][idx.sup])-1
     
     # cutoff for median
     sup.med <- median(sup.pred.sorted)
     idx.sup.med <- which.min(sup.pred.sorted<sup.med) # index at midrange
     
     # cutoff determination for PP
     idx.sup.pp50 <- which.min(sup.pred.sorted<0.5) # index at pp50 cutoff
     idx.sup.pp80 <- which.min(sup.pred.sorted<0.8) # index at pp80 cutoff
     idx.sup.pp95 <- which.min(sup.pred.sorted<0.95) # index at pp95 cutoff
     
     # number of positive predicted  
     sup.pos.med[i,j] <- sum(sup.lab.sorted[idx.sup.med:length(sup.lab.sorted)])
    
     sup.pos.pp50[i,j] <- sum(sup.lab.sorted[idx.sup.pp50:length(sup.lab.sorted)])
     sup.pos.pp80[i,j] <- sum(sup.lab.sorted[idx.sup.pp80:length(sup.lab.sorted)])
     sup.pos.pp95[i,j] <- sum(sup.lab.sorted[idx.sup.pp95:length(sup.lab.sorted)])
     
     # number of predicted positive saved by contamination rate
     if(g==1){
       sup.pos.med.0 <- sup.pos.med
       
       sup.pos.pp50.0 <- sup.pos.pp50
       sup.pos.pp80.0 <- sup.pos.pp80
       sup.pos.pp95.0 <- sup.pos.pp95
     }
     if(g==3){
       sup.pos.med.20 <- sup.pos.med
       
       sup.pos.pp50.20 <- sup.pos.pp50
       sup.pos.pp80.20 <- sup.pos.pp80
       sup.pos.pp95.20 <- sup.pos.pp95
     }
     if(g==6){
       sup.pos.med.50 <- sup.pos.med
       
       sup.pos.pp50.50 <- sup.pos.pp50
       sup.pos.pp80.50 <- sup.pos.pp80
       sup.pos.pp95.50 <- sup.pos.pp95
     }

     # precision
     sup.prec.med[i,j] <- sum(sup.lab.sorted[idx.sup.med:length(sup.lab.sorted)])/length(idx.sup.med:length(sup.lab.sorted))
     
     sup.prec.pp50[i,j] <- sum(sup.lab.sorted[idx.sup.pp50:length(sup.lab.sorted)])/length(idx.sup.pp50:length(sup.lab.sorted))
     sup.prec.pp80[i,j] <- sum(sup.lab.sorted[idx.sup.pp80:length(sup.lab.sorted)])/length(idx.sup.pp80:length(sup.lab.sorted))
     sup.prec.pp95[i,j] <- sum(sup.lab.sorted[idx.sup.pp95:length(sup.lab.sorted)])/length(idx.sup.pp95:length(sup.lab.sorted))
     
     # recall
     sup.rec.med[i,j] <- sum(sup.lab.sorted[idx.sup.med:length(sup.lab.sorted)])/sum(sup.lab.sorted)
     
     sup.rec.pp50[i,j] <- sum(sup.lab.sorted[idx.sup.pp50:length(sup.lab.sorted)])/sum(sup.lab.sorted)
     sup.rec.pp80[i,j] <- sum(sup.lab.sorted[idx.sup.pp80:length(sup.lab.sorted)])/sum(sup.lab.sorted)
     sup.rec.pp95[i,j] <- sum(sup.lab.sorted[idx.sup.pp95:length(sup.lab.sorted)])/sum(sup.lab.sorted)
     
     # png(paste0("results/plots/Histograms/LASSO_Histogram_", p[i]*100,"_Contamination_", contam.perc[g]*100,".png"),height = 8, width = 8, res = 300, units = 'in')
     # plot(density(sup.pred), main = paste0("Kernal Density of LASSO Scaled Scores \nTraining Size ", p[i]*100, " % Contamination ", contam.perc[g]*100,"%"),xlab = "Probability",xlim=c(0,1))
     # dev.off()
     }
     else{
       # precision
       sup.prec.med[i,j] <- 0
       
       sup.prec.pp50[i,j] <- 0
       sup.prec.pp80[i,j] <- 0
       sup.prec.pp95[i,j] <- 0
       
       # recall
       sup.rec.med[i,j] <- 0
       
       sup.rec.pp50[i,j] <- 0
       sup.rec.pp80[i,j] <- 0
       sup.rec.pp95[i,j] <- 0
       
       
     }
     
# Semi performance cutoffs
     if(g<2){

       semi.pred1 <- prediction(all.preds.semi[[i]],all.labels.semi[[i]]) # predicted response vectors
       # semi.pred <- slot(semi.pred1,'predictions')[[1]]/max(slot(semi.pred1,'predictions')[[1]])
       
       semi.pred <- (slot(semi.pred1,'predictions')[[1]]-min(slot(semi.pred1,'predictions')[[1]]))/(max(slot(semi.pred1,'predictions')[[1]])-min(slot(semi.pred1,'predictions')[[1]]))
       summary(semi.pred)
       idx.semi <- order(semi.pred)
       semi.pred.sorted <- semi.pred[idx.semi]
       semi.lab.sorted <- as.numeric(slot(semi.pred1,'labels')[[1]][idx.semi])-1
       
       # cutoff for median
       semi.med <- median(semi.pred.sorted)
       idx.semi.med <- which.min(semi.pred.sorted<semi.med) # index at midrange
       
       # cutoff determination for PP
       idx.semi.pp50 <- which.min(semi.pred.sorted<0.5) # index at pp50 cutoff
       idx.semi.pp80 <- which.min(semi.pred.sorted<0.8) # index at pp80 cutoff
       idx.semi.pp95 <- which.min(semi.pred.sorted<0.95) # index at pp95 cutoff
       
       # number of positive predicted  
       semi.pos.med[i,j] <- sum(semi.lab.sorted[idx.semi.med:length(semi.lab.sorted)])
       
       semi.pos.pp50[i,j] <- sum(semi.lab.sorted[idx.semi.pp50:length(semi.lab.sorted)])
       semi.pos.pp80[i,j] <- sum(semi.lab.sorted[idx.semi.pp80:length(semi.lab.sorted)])
       semi.pos.pp95[i,j] <- sum(semi.lab.sorted[idx.semi.pp95:length(semi.lab.sorted)])
       
       # number of predicted positive saved by contamination rate
       if(g==1){
         semi.pos.med.0 <- semi.pos.med
         
         semi.pos.pp50.0 <- semi.pos.pp50
         semi.pos.pp80.0 <- semi.pos.pp80
         semi.pos.pp95.0 <- semi.pos.pp95
       }
       if(g==3){
         semi.pos.med.20 <- semi.pos.med
         
         semi.pos.pp50.20 <- semi.pos.pp50
         semi.pos.pp80.20 <- semi.pos.pp80
         semi.pos.pp95.20 <- semi.pos.pp95
       }
       if(g==6){
         semi.pos.med.50 <- semi.pos.med
         
         semi.pos.pp50.50 <- semi.pos.pp50
         semi.pos.pp80.50 <- semi.pos.pp80
         semi.pos.pp95.50 <- semi.pos.pp95
       }
       
       # precision
       semi.prec.med[i,j] <- sum(semi.lab.sorted[idx.semi.med:length(semi.lab.sorted)])/length(idx.semi.med:length(semi.lab.sorted))
       
       semi.prec.pp50[i,j] <- sum(semi.lab.sorted[idx.semi.pp50:length(semi.lab.sorted)])/length(idx.semi.pp50:length(semi.lab.sorted))
       semi.prec.pp80[i,j] <- sum(semi.lab.sorted[idx.semi.pp80:length(semi.lab.sorted)])/length(idx.semi.pp80:length(semi.lab.sorted))
       semi.prec.pp95[i,j] <- sum(semi.lab.sorted[idx.semi.pp95:length(semi.lab.sorted)])/length(idx.semi.pp95:length(semi.lab.sorted))
       
       # recall
       semi.rec.med[i,j] <- sum(semi.lab.sorted[idx.semi.med:length(semi.lab.sorted)])/sum(semi.lab.sorted)
       
       semi.rec.pp50[i,j] <- sum(semi.lab.sorted[idx.semi.pp50:length(semi.lab.sorted)])/sum(semi.lab.sorted)
       semi.rec.pp80[i,j] <- sum(semi.lab.sorted[idx.semi.pp80:length(semi.lab.sorted)])/sum(semi.lab.sorted)
       semi.rec.pp95[i,j] <- sum(semi.lab.sorted[idx.semi.pp95:length(semi.lab.sorted)])/sum(semi.lab.sorted)
       
    # png(paste0("results/plots/Histograms/Semi_Histogram_", p[i]*100,"_Contamination_", contam.perc[g]*100,".png"),height = 8, width = 8, res = 300, units = 'in')
    # plot(density(semi.pred),main = paste0("Kernal Density of Semi-Supervised Predicted Probabilities\nTraining Size ", p[i]*100, " % Contamination ", contam.perc[g]*100,"%"),xlab = "Probability",xlim=c(0,1))
    # dev.off()
     }
    }
    }  
    
    if(g==c(1,3,6) & j==1){
      for(i in 1:l){

      lab <- factor(all.labels.sup[[i]], labels = c("Positive","Negative"),levels = c(1,0))
      score <- as.vector(all.preds.sup.p[[i]])
      pos.lab.idx <- which(lab=='Positive')
      bw <- 0.05
      png(paste0("results/plots/Histograms/Lasso_Histogram_", p[i]*100,"_Training_", contam.perc[g]*100,"_Contamination.png"),height = 8, width = 8, res = 300, units = 'in')
      plot(density(score[pos.lab.idx],bw=bw), xlab='Predicted Score', ylab='Density',ylim=c(0,5),cex.axis=1.5,font=2,main="",cex.lab=1.5,font.lab=2)
      lines(density(score[-pos.lab.idx],bw=bw),lty=2)
      legend('top', levels(lab),lty = c(1,2), cex = 1.5,text.font=2)
      dev.off()
      
      lab1 <- factor(all.labels.semi[[i]], labels = c("Positive","Negative"),levels = c(1,0))
      score1 <- as.vector(all.preds.semi[[i]])
      pos.lab.idx <- which(lab1=='Positive')
      bw <- 0.05
      png(paste0("results/plots/Histograms/Semi_Histogram_", p[i]*100,"_Training_", contam.perc[g]*100,"_Contamination.png"),height = 8, width = 8, res = 300, units = 'in')
      plot(density(score1[pos.lab.idx],bw=bw), xlab='Predicted Score', ylab='Density',ylim=c(0,4),cex.axis=1.5,font=2,main="",cex.lab=1.5,font.lab=2)
      lines(density(score1[-pos.lab.idx],bw=bw),lty=2)
      legend('top', levels(lab),lty = c(1,2), cex = 1.5,text.font=2)
      dev.off()
                   }
                }
    
    
    ROC.svm <- multiroc(all.preds.svm,all.labels.svm) # collection of ROC calculations from prediction and label lists
    ROC.J48 <- multiroc(all.preds.J48,all.labels.J48) # collection of ROC calculations from prediction and label lists
    ROC.sup <- multiroc(all.preds.sup.p,all.labels.sup) # collection of ROC calculations from prediction and label lists
  if(g<2){ROC.semi <- multiroc(all.preds.semi,all.labels.semi)
         perf.semi <- all.perfs.semi}
   
    percent.contam = list() # list of contamination percentages
    
    for(m in 1:length(number.pos)){
      percent.contam$new <- number.contam[m]/(inTrain.neg.only.length[m]+number.contam[m]) # contamination percentages determined from number of contaminates/total negative labels
      names(percent.contam)[length(percent.contam)] = p[m] # assigns names to list entries 
    }
    
    # Find AUC and collect contamination percents
    # Each row represents a training percent
    # Each column represents an iteration
    for(k in 1:l){
      AUC.svm[k,j] <- ROC.svm[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
      label.svm[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
      AUC.J48[k,j] <- ROC.J48[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
      label.tree[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
      AUC.sup[k,j] <- ROC.sup[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
      label.lasso[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
      
      semi.mean <-apply(AUC.semi,1,mean)
      semi.median <-apply(AUC.semi,1,median)
      semi.var <-apply(AUC.semi,1,var)
      semi.mad <- apply(AUC.semi,1,mad)
      
      svm.mean[[g]] <-apply(AUC.svm,1,mean)
      svm.median[[g]] <-apply(AUC.svm,1,median)
      svm.var[[g]] <-apply(AUC.svm,1,var)
      svm.mad[[g]] <-apply(AUC.svm,1,mad)
      forest.mean[[g]] <-apply(AUC.J48,1,mean)
      forest.median[[g]] <-apply(AUC.J48,1,median)
      forest.var[[g]] <-apply(AUC.J48,1,var)
      forest.mad[[g]] <-apply(AUC.J48,1,mad)
      sup.mean[[g]] <-apply(AUC.sup,1,mean)
      sup.median[[g]] <-apply(AUC.sup,1,median)
      sup.var[[g]] <-apply(AUC.sup,1,var)
      sup.mad[[g]] <-apply(AUC.sup,1,mad)
      
      if(g<2){AUC.semi[k,j] <- ROC.semi[[k]]$auc}
    }
  }  
  
  svm.contam <- apply(label.svm,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  tree.contam <- apply(label.tree,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  sup.contam <- apply(label.lasso,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  
  AUC2 <- t(AUC.svm)
  AUC1 <- t(AUC.semi)
  AUC3 <- t(AUC.sup)
  AUC4 <- t(AUC.J48)
  AUC <- matrix(0,ncol=4*l,nrow=n1)
  
  for(f in 1:ncol(AUC1)){
    AUC[,4*f-3] <- AUC1[,f]
    AUC[,4*f-2] <- AUC2[,f]
    AUC[,4*f-1] <- AUC3[,f]
    AUC[,4*f]   <- AUC4[,f]
  }
 

  AUC.list$new <- as.vector(AUC)
  names(AUC.list)[length(AUC.list)] = paste0(contam.perc[g]*100,'% contamination')

# for(g in 1:length(contam.perc)){  # Use for rerunning plots only
  n <- rep(p,each=4*n1)
  n <- factor(n)
  method <- rep(c(rep('semi',n1),
                  rep('svm',n1),
                  rep('sup',n1),
                  rep('tree',n1)),l)
  
  tog <- data.frame(AUC=AUC.list[[g]],n=n,method=method)
  box.at <- seq(1:(l*5))
  remove <- seq(5,length(box.at),5)
  at <- box.at[!(box.at%in%remove)]
  
  png(paste0('results/plots/6D_allFeatures_unbalanced_',contam.perc[g]*100,'.png'),
      height = 8, width = 8, res = 300, units = 'in')
  boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','dodgerblue2','cadetblue3','lightblue1'),xlab='Training Set (%)', ylab='AUC',xaxt="n",
          main=paste0('Supervised vs Semi-supervised (All 17 Features)\nAUC for various training % with ',n1,' iterations\nTest: unbalanced (',contam.perc[g]*100,'% Contamination)'),ylim=c(0.38,0.8),cex.axis=1.5,cex.lab=1.4)
  axis(side=1,at=seq(2.5,92.5,5),labels=seq(1,10,0.5),las=2,cex.axis=1.5)
  legend("topleft",c("Semi-Supervised","Supervised-Lasso","Supervised-SVM","Supervised-SVM"),col=c('red','dodgerblue2','cadetblue3','lightblue1'),pch=15)
  dev.off()

  
  # For paper without title
  png(paste0('results/plots/6D_allFeatures_unbalanced_',contam.perc[g]*100,'_paper.png'),
      height = 8, width = 8, res = 300, units = 'in')
  boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','dodgerblue2','cadetblue3','lightblue1'),xlab='Training Set (%)', ylab='AUC',xaxt="n",ylim=c(0.38,0.8),cex.axis=1.5,cex.lab=1.4)
  axis(side=1,at=seq(2.5,92.5,5),labels=seq(1,10,0.5),las=2,cex.axis=1.5)
  legend("topleft",c("Semi-Supervised","Supervised-Lasso","Supervised-SVM","Supervised-Forest"),col=c('red','dodgerblue2','cadetblue3','lightblue1'),pch=15)
  dev.off()

 # } # Use for running AUC plots only

if(g==1|g==3|g==6){
# precision tp/(tp+fp) = PPV
  # cutoff at median
  sup.precision.med <- apply(sup.prec.med,1,function(x) mean(x,na.rm = T))
  names(sup.precision.med) <- p
  semi.precision.med <- apply(semi.prec.med,1,function(x) mean(x,na.rm = T))
  names(semi.precision.med) <- p
  prec0.med <- rbind(sup.precision.med,semi.precision.med)
  
  # cutoff at midrange
  sup.precision.midrange <- apply(sup.prec.midrange,1,function(x) mean(x,na.rm = T))
  names(sup.precision.midrange) <- p
  semi.precision.midrange <- apply(semi.prec.midrange,1,function(x) mean(x,na.rm = T))
  names(semi.precision.midrange) <- p
  prec0.midrange <- rbind(sup.precision.midrange,semi.precision.midrange)
  
  # cutoff at pp50
  sup.precision.pp50 <- apply(sup.prec.pp50,1,function(x) mean(x,na.rm = T))
  names(sup.precision.pp50) <- p
  semi.precision.pp50 <- apply(semi.prec.pp50,1,function(x) mean(x,na.rm = T))
  names(semi.precision.pp50) <- p
  prec0.pp50 <- rbind(sup.precision.pp50,semi.precision.pp50)
  
  # cutoff at pp80
  sup.precision.pp80 <- apply(sup.prec.pp80,1,function(x) mean(x,na.rm = T))
  names(sup.precision.pp80) <- p
  semi.precision.pp80 <- apply(semi.prec.pp80,1,function(x) mean(x,na.rm = T))
  names(semi.precision.pp80) <- p
  prec0.pp80 <- rbind(sup.precision.pp80,semi.precision.pp80)
  
  # cutoff at pp95
  sup.precision.pp95 <- apply(sup.prec.pp95,1,function(x) mean(x,na.rm = T))
  names(sup.precision.pp90) <- p
  semi.precision.pp95 <- apply(semi.prec.pp95,1,function(x) mean(x,na.rm = T))
  names(semi.precision.pp95) <- p
  prec0.pp95 <- rbind(sup.precision.pp95,semi.precision.pp95)
  
# recall tp/(tp+fn) 
  # cutoff at median
  sup.recall.med <- apply(sup.rec.med,1,function(x) mean(x,na.rm = T))
  names(sup.recall.med) <- p
  semi.recall.med <- apply(semi.rec.med,1,function(x) mean(x,na.rm = T))
  names(semi.recall.med) <- p
  rec0.med <- rbind(sup.recall.med,semi.recall.med)
  
  # cutoff at midrange
  sup.recall.midrange <- apply(sup.rec.midrange,1,function(x) mean(x,na.rm = T))
  names(sup.recall.midrange) <- p
  semi.recall.midrange <- apply(semi.rec.midrange,1,function(x) mean(x,na.rm = T))
  names(semi.recall.midrange) <- p
  rec0.midrange <- rbind(sup.recall.midrange,semi.recall.midrange)
  
  # cutoff at pp50
  sup.recall.pp50 <- apply(sup.rec.pp50,1,function(x) mean(x,na.rm = T))
  names(sup.recall.pp50) <- p
  semi.recall.pp50 <- apply(semi.rec.pp50,1,function(x) mean(x,na.rm = T))
  names(semi.recall.pp50) <- p
  rec0.pp50 <- rbind(sup.recall.pp50,semi.recall.pp50)
  
  # cutoff at pp80
  sup.recall.pp80 <- apply(sup.rec.pp80,1,function(x) mean(x,na.rm = T))
  names(sup.recall.pp80) <- p
  semi.recall.pp80 <- apply(semi.rec.pp80,1,function(x) mean(x,na.rm = T))
  names(semi.recall.pp80) <- p
  rec0.pp80 <- rbind(sup.recall.pp80,semi.recall.pp80)
  
  # cutoff at pp90
  sup.recall.pp95 <- apply(sup.rec.pp95,1,function(x) mean(x,na.rm = T))
  names(sup.recall.pp95) <- p
  semi.recall.pp95 <- apply(semi.rec.pp95,1,function(x) mean(x,na.rm = T))
  names(semi.recall.pp95) <- p
  rec0.pp95 <- rbind(sup.recall.pp95,semi.recall.pp95)
  }
  
if(g==1){  
  # Table for performance measures 
  # cutoff at median
  summary0.med <- rbind(prec0.med[1,],rec0.med[1,],f0.med[1,],prec0.med[2,],rec0.med[2,],f0.med[2,])
  rownames(summary0.med) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
  colnames(summary0.med) <- seq(1,10,0.5)
  summary0.med[3,] <- apply(summary0.med,2,function(x) 2/((1/x[1])+(1/x[2])))
  summary0.med[6,] <- apply(summary0.med,2,function(x) 2/((1/x[4])+(1/x[5])))
  write.csv(summary0.med,'results/tables/perf0.med.csv')
  
  # cutoff at midrange 
  summary0.midrange <- rbind(prec0.midrange[1,],rec0.midrange[1,],f0.midrange[1,],prec0.midrange[2,],rec0.midrange[2,],f0.midrange[2,])
  rownames(summary0.midrange) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
  colnames(summary0.midrange) <- seq(1,10,0.5)
  summary0.midrange[3,] <- apply(summary0.midrange,2,function(x) 2/((1/x[1])+(1/x[2])))
  summary0.midrange[6,] <- apply(summary0.midrange,2,function(x) 2/((1/x[4])+(1/x[5])))
  write.csv(summary0.midrange,'results/tables/perf0.midrange.csv')
  
  # cutoff at pp50
  summary0.pp50 <- rbind(prec0.pp50[1,],rec0.pp50[1,],f0.pp50[1,],prec0.pp50[2,],rec0.pp50[2,],f0.pp50[2,])
  rownames(summary0.pp50) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
  colnames(summary0.pp50) <- seq(1,10,0.5)
  summary0.pp50[3,] <- apply(summary0.pp50,2,function(x) 2/((1/x[1])+(1/x[2])))
  summary0.pp50[6,] <- apply(summary0.pp50,2,function(x) 2/((1/x[4])+(1/x[5])))
  write.csv(summary0.pp50,'results/tables/perf0.pp50.csv')
  
  # cutoff at pp80
  summary0.pp80 <- rbind(prec0.pp80[1,],rec0.pp80[1,],f0.pp80[1,],prec0.pp80[2,],rec0.pp80[2,],f0.pp80[2,])
  rownames(summary0.pp80) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
  colnames(summary0.pp80) <- seq(1,10,0.5)
  summary0.pp80[3,] <- apply(summary0.pp80,2,function(x) 2/((1/x[1])+(1/x[2])))
  summary0.pp80[6,] <- apply(summary0.pp80,2,function(x) 2/((1/x[4])+(1/x[5])))
  write.csv(summary0.pp80,'results/tables/perf0.pp80.csv')
  
  # cutoff at pp90
  summary0.pp95 <- rbind(prec0.pp95[1,],rec0.pp95[1,],f0.pp95[1,],prec0.pp95[2,],rec0.pp95[2,],f0.pp95[2,])
  rownames(summary0.pp95) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
  colnames(summary0.pp95) <- seq(1,10,0.5)
  summary0.pp95[3,] <- apply(summary0.pp95,2,function(x) 2/((1/x[1])+(1/x[2])))
  summary0.pp95[6,] <- apply(summary0.pp95,2,function(x) 2/((1/x[4])+(1/x[5])))
  write.csv(summary0.pp95,'results/tables/perf0.pp95.csv')
}

  if(g==3){

    # cutoff at median
    summary20.med <- rbind(prec0.med[1,],rec0.med[1,],f0.med[1,],prec0.med[2,],rec0.med[2,],f0.med[2,])
    rownames(summary20.med) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary20.med) <- seq(1,10,0.5)
    summary20.med[3,] <- apply(summary20.med,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary20.med[6,] <- apply(summary20.med,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary20.med,'results/tables/perf20.med.csv')
    
    # cutoff at midrange 
    summary20.midrange <- rbind(prec0.midrange[1,],rec0.midrange[1,],f0.midrange[1,],prec0.midrange[2,],rec0.midrange[2,],f0.midrange[2,])
    rownames(summary20.midrange) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary20.midrange) <- seq(1,10,0.5)
    summary20.midrange[3,] <- apply(summary20.midrange,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary20.midrange[6,] <- apply(summary20.midrange,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary20.midrange,'results/tables/perf20.midrange.csv')
    
    # cutoff at pp50
    summary20.pp50 <- rbind(prec0.pp50[1,],rec0.pp50[1,],f0.pp50[1,],prec0.pp50[2,],rec0.pp50[2,],f0.pp50[2,])
    rownames(summary20.pp50) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary20.pp50) <- seq(1,10,0.5)
    summary20.pp50[3,] <- apply(summary20.pp50,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary20.pp50[6,] <- apply(summary20.pp50,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary20.pp50,'results/tables/perf20.pp50.csv')
    
    # cutoff at pp80
    summary20.pp80 <- rbind(prec0.pp80[1,],rec0.pp80[1,],f0.pp80[1,],prec0.pp80[2,],rec0.pp80[2,],f0.pp80[2,])
    rownames(summary20.pp80) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary20.pp80) <- seq(1,10,0.5)
    summary20.pp80[3,] <- apply(summary20.pp80,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary20.pp80[6,] <- apply(summary20.pp80,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary20.pp80,'results/tables/perf20.pp80.csv')
    
    # cutoff at pp95
    summary20.pp95 <- rbind(prec0.pp95[1,],rec0.pp95[1,],f0.pp95[1,],prec0.pp95[2,],rec0.pp95[2,],f0.pp95[2,])
    rownames(summary20.pp95) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary20.pp95) <- seq(1,10,0.5)
    summary20.pp95[3,] <- apply(summary20.pp95,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary20.pp95[6,] <- apply(summary20.pp95,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary20.pp95,'results/tables/perf20.pp95.csv')
    }
  
  if(g==6){

    # cutoff at median
    summary50.med <- rbind(prec0.med[1,],rec0.med[1,],f0.med[1,],prec0.med[2,],rec0.med[2,],f0.med[2,])
    rownames(summary50.med) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary50.med) <- seq(1,10,0.5)
    summary50.med[3,] <- apply(summary50.med,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary50.med[6,] <- apply(summary50.med,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary50.med,'results/tables/perf50.med.csv')
    
    # cutoff at midrange 
    summary50.midrange <- rbind(prec0.midrange[1,],rec0.midrange[1,],f0.midrange[1,],prec0.midrange[2,],rec0.midrange[2,],f0.midrange[2,])
    rownames(summary50.midrange) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary50.midrange) <- seq(1,10,0.5)
    summary50.midrange[3,] <- apply(summary50.midrange,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary50.midrange[6,] <- apply(summary50.midrange,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary50.midrange,'results/tables/perf50.midrange.csv')
    
    # cutoff at pp50
    summary50.pp50 <- rbind(prec0.pp50[1,],rec0.pp50[1,],f0.pp50[1,],prec0.pp50[2,],rec0.pp50[2,],f0.pp50[2,])
    rownames(summary50.pp50) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary50.pp50) <- seq(1,10,0.5)
    summary50.pp50[3,] <- apply(summary50.pp50,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary50.pp50[6,] <- apply(summary50.pp50,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary50.pp50,'results/tables/perf50.pp50.csv')
    
    # cutoff at pp80
    summary50.pp80 <- rbind(prec0.pp80[1,],rec0.pp80[1,],f0.pp80[1,],prec0.pp80[2,],rec0.pp80[2,],f0.pp80[2,])
    rownames(summary50.pp80) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary50.pp80) <- seq(1,10,0.5)
    summary50.pp80[3,] <- apply(summary50.pp80,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary50.pp80[6,] <- apply(summary50.pp80,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary50.pp80,'results/tables/perf50.pp80.csv')
    
    # cutoff at pp90
    summary50.pp95 <- rbind(prec0.pp95[1,],rec0.pp95[1,],f0.pp95[1,],prec0.pp95[2,],rec0.pp95[2,],f0.pp95[2,])
    rownames(summary50.pp95) <- c('precision sup','recall sup','f-measure sup','precision Semi','recall Semi','f-measure Semi')
    colnames(summary50.pp95) <- seq(1,10,0.5)
    summary50.pp95[3,] <- apply(summary50.pp95,2,function(x) 2/((1/x[1])+(1/x[2])))
    summary50.pp95[6,] <- apply(summary50.pp95,2,function(x) 2/((1/x[4])+(1/x[5])))
    write.csv(summary50.pp95,'results/tables/perf50.pp95.csv')
  }
  }
  
method.means <- cbind(semi.mean,svm.mean[[1]],forest.mean[[1]],sup.mean[[1]],svm.mean[[2]],forest.mean[[2]],sup.mean[[2]],svm.mean[[3]],forest.mean[[3]],sup.mean[[3]])
method.median <- cbind(semi.median,svm.median[[1]],forest.median[[1]],sup.median[[1]],svm.median[[2]],forest.median[[2]],sup.median[[2]],svm.median[[3]],forest.median[[3]],sup.median[[3]])
method.var <- cbind(semi.var,svm.var[[1]],forest.var[[1]],sup.var[[1]],svm.var[[2]],forest.var[[2]],sup.var[[2]],svm.var[[3]],forest.var[[3]],sup.var[[3]])
method.mad <- cbind(semi.mad,svm.mad[[1]],forest.mad[[1]],sup.mad[[1]],svm.mad[[2]],forest.mad[[2]],sup.mad[[2]],svm.mad[[3]],forest.mad[[3]],sup.mad[[3]])

colnames(method.means) <- c('Semi-sup','SVM 0 %','Forest 0 %','Lasso 0 %','SVM 20 %','Forest 20 %','Lasso 20 %','SVM 50 %','Forest 50 %','Lasso 50 %')
colnames(method.median) <- c('Semi-sup','SVM 0 %','Forest 0 %','Lasso 0 %','SVM 20 %','Forest 20 %','Lasso 20 %','SVM 50 %','Forest 50 %','Lasso 50 %')
colnames(method.var) <- c('Semi-sup','SVM 0 %','Forest 0 %','Lasso 0 %','SVM 20 %','Forest 20 %','Lasso 20 %','SVM 50 %','Forest 50 %','Lasso 50 %')
colnames(method.var) <- c('Semi-sup','SVM 0 %','Forest 0 %','Lasso 0 %','SVM 20 %','Forest 20 %','Lasso 20 %','SVM 50 %','Forest 50 %','Lasso 50 %')
write.csv(method.means,"results/tables/boxplot_means.csv")
write.csv(method.median,"results/tables/boxplot_median.csv")
write.csv(method.var,"results/tables/boxplot_var.csv")
write.csv(method.mad,"results/tables/boxplot_mad.csv")

perc <- rep(as.numeric(rownames(method.means)),4)*100
mean0 <- c(method.means[,1],method.means[,2],method.means[,3],method.means[,4])
mean20 <- c(method.means[,1],method.means[,5],method.means[,6],method.means[,7])
mean50 <- c(method.means[,1],method.means[,8],method.means[,9],method.means[,10])
variance0 <- c(method.var[,1],method.var[,2],method.var[,3],method.var[,4])
variance20 <- c(method.var[,1],method.var[,5],method.var[,6],method.var[,7])
variance50 <- c(method.var[,1],method.var[,8],method.var[,9],method.var[,10])
median0 <- c(method.median[,1],method.median[,2],method.median[,3],method.median[,4])
median20 <- c(method.median[,1],method.median[,5],method.median[,6],method.median[,7])
median50 <- c(method.median[,1],method.median[,8],method.median[,9],method.median[,10])
mad0 <- c(method.mad[,1],method.mad[,2],method.mad[,3],method.mad[,4])
mad20 <- c(method.mad[,1],method.mad[,5],method.mad[,6],method.mad[,7])
mad50 <- c(method.mad[,1],method.mad[,8],method.mad[,9],method.mad[,10])
cv0 <- sqrt(variance0)/mean0
cv20 <- sqrt(variance20)/mean20
cv50 <- sqrt(variance50)/mean50
grp <- rep(c('Semi-sup','SVM','Forest','Lasso'),each=length(method.var[,1]))
df0.mean <- data.frame(perc,mean0,Method=grp)
df20.mean <- data.frame(perc,mean20,Method=grp)
df50.mean <- data.frame(perc,mean50,Method=grp)
df0.median <- data.frame(perc,median0,Method=grp)
df20.median <- data.frame(perc,median20,Method=grp)
df50.median <- data.frame(perc,median50,Method=grp)
df0.var <- data.frame(perc,variance0,Method=grp)
df20.var <- data.frame(perc,variance20,Method=grp)
df50.var <- data.frame(perc,variance50,Method=grp)
df0.mad <- data.frame(perc,mad0,Method=grp)
df20.mad <- data.frame(perc,mad20,Method=grp)
df50.mad <- data.frame(perc,mad50,Method=grp)
df0.cv.med <- data.frame(perc,cv=mad0/median0,Method=grp)
df20.cv.med <- data.frame(perc,cv=mad20/median20,Method=grp)
df50.cv.med <- data.frame(perc,cv=mad50/median50,Method=grp)
df0.cv.mean <- data.frame(perc,cv=cv0,Method=grp)
df20.cv.mean <- data.frame(perc,cv=cv20,Method=grp)
df50.cv.mean <- data.frame(perc,cv=cv50,Method=grp)

levels(df0.mean$Method) <- levels(df0.var$Method) <- levels(df20.mean$Method) <- levels(df20.var$Method) <-
  levels(df50.mean$Method) <- levels(df50.var$Method) <- levels(df0.cv$Method) <- levels(df20.cv$Method) <-
  levels(df50.cv.med$Method) <- levels(df50.cv.mean$Method) <- c('Forest','Lasso','Semi-sup','SVM')

### Plots for AUC means and variances
p1 <-ggplot(df0.var, aes(x=perc,y=variance0,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Variance',x='Training Set (%)')) + 
  ylim(0,0.0065) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())


p2 <- ggplot(df20.var, aes(x=perc,y=variance20,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Variance',x='Training Set (%)')) + 
  ylim(0,0.0065) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p3 <- ggplot(df50.var, aes(x=perc,y=variance50,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Variance',x='Training Set (%)')) + 
  ylim(0,0.0065) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())


png('results/plots/AUCvariances.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p1,p2,p3)
dev.off()

p4 <-ggplot(df0.mean, aes(x=perc,y=mean0,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Mean',x='Training Set (%)')) +
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p5 <- ggplot(df20.mean, aes(x=perc,y=mean20,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Mean',x='Training Set (%)')) + 
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p6 <- ggplot(df50.mean, aes(x=perc,y=mean50,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Mean',x='Training Set (%)')) + 
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

png('results/plots/AUCmeans.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4,p5,p6)
dev.off()

p4.cv.mean <-ggplot(df0.cv.mean, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC SD/Mean',x='Training Set (%)')) +
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = c(0.65,0.70)) + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p5.cv.mean <- ggplot(df20.cv.mean, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC SD/Mean',x='Training Set (%)')) + 
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p6.cv.mean <- ggplot(df50.cv.mean, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC SD/Mean',x='Training Set (%)')) + 
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())


png('results/plots/AUCcvmean.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4.cv.mean,p5.cv.mean,p6.cv.mean)
dev.off()

# mad, median, CV
p4.med <-ggplot(df0.median, aes(x=perc,y=median0,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC median',x='Training Set (%)')) +
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p5.med <- ggplot(df20.median, aes(x=perc,y=median20,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC median',x='Training Set (%)')) + 
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p6.med <- ggplot(df50.median, aes(x=perc,y=median50,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC Median',x='Training Set (%)')) + 
  ylim(0.45,0.7) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())

png('results/plots/AUCmed.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4.med,p5.med,p6.med)
dev.off()

p4.mad <-ggplot(df0.mad, aes(x=perc,y=mad0,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD',x='Training Set (%)')) +
  ylim(0,0.15) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p5.mad <- ggplot(df20.mad, aes(x=perc,y=mad20,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD',x='Training Set (%)')) + 
  ylim(0,0.15) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12)) + theme(legend.key=element_blank())

p6.mad <- ggplot(df50.mad, aes(x=perc,y=mad50,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD',x='Training Set (%)')) + 
  ylim(0,0.15) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12)) + theme(legend.key=element_blank())

png('results/plots/AUCmad.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4.mad,p5.mad,p6.mad)
dev.off()

p4.cv.med <-ggplot(df0.cv.med, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD/Median',x='Training Set (%)')) +
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = c(0.65,0.70)) + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) + theme(legend.key=element_blank())
p5.cv.med <- ggplot(df20.cv.med, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD/Median',x='Training Set (%)')) + 
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12))
p6.cv.med <- ggplot(df50.cv.med, aes(x=perc,y=cv,group=Method,colour=Method)) + geom_line(size=1) + 
  labs(list(y='AUC MAD/Median',x='Training Set (%)')) + 
  ylim(0,0.25) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('lightblue1','dodgerblue2','red','cadetblue3')) +
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12))

png('results/plots/AUCcvmed.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4.cv.med,p5.cv.med,p6.cv.med)
dev.off()

# Transpose mean, variance, and CV plots

png('results/plots/perf0.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p4,p1,p4.cv)
dev.off()

png('results/plots/perf20.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p5,p2,p5.cv)
dev.off()

png('results/plots/perf50.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p6,p3,p6.cv)
dev.off()

### Plots for precision, recall, and f
# Median cutoff
performance.vector.0 <- c(summary0.med[1,],summary0.med[2,],summary0.med[3,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
performance.vector.20 <- c(summary20.med[1,],summary20.med[2,],summary20.med[3,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
performance.vector.50 <- c(summary50.med[1,],summary50.med[2,],summary50.med[3,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
perc.1 <- rep(as.numeric(rownames(method.means)),6)*100
measure <- rep(c(rep('Precision',19),rep('Recall',19),rep('f-measure',19)),2)
method.1 <- c(rep('Lasso',57),rep('Semi',57))
group <- interaction(measure,method.1)
median.0 <- data.frame(percent=perc.1,performance=performance.vector.0,measure=measure,method=method.1,group=group)
median.20 <- data.frame(percent=perc.1,performance=performance.vector.20,measure=measure,method=method.1,group=group)
median.50 <- data.frame(percent=perc.1,performance=performance.vector.50,measure=measure,method=method.1,group=group)

p7 <- ggplot(median.0, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=16,face = 'bold'),
      axis.title=element_text(size=16,face = 'bold'))

p8 <- ggplot(median.20, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

p9 <- ggplot(median.50, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

png('results/plots/measures_median.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7,p8,p9)
dev.off()

# pp50 cutoff
performance.vector.0 <- c(summary0.pp50[1,],summary0.pp50[2,],summary0.pp50[3,],summary0.pp50[4,],summary0.pp50[5,],summary0.pp50[6,])
performance.vector.20 <- c(summary20.pp50[1,],summary20.pp50[2,],summary20.pp50[3,],summary0.pp50[4,],summary0.pp50[5,],summary0.pp50[6,])
performance.vector.50 <- c(summary50.pp50[1,],summary50.pp50[2,],summary50.pp50[3,],summary0.pp50[4,],summary0.pp50[5,],summary0.pp50[6,])
perc.1 <- rep(as.numeric(rownames(method.means)),6)*100
measure <- rep(c(rep('Precision',19),rep('Recall',19),rep('f-measure',19)),2)
method.1 <- c(rep('Lasso',57),rep('Semi',57))
group <- interaction(measure,method.1)
median.0 <- data.frame(percent=perc.1,performance=performance.vector.0,measure=measure,method=method.1,group=group)
median.20 <- data.frame(percent=perc.1,performance=performance.vector.20,measure=measure,method=method.1,group=group)
median.50 <- data.frame(percent=perc.1,performance=performance.vector.50,measure=measure,method=method.1,group=group)

p7.pp50 <- ggplot(median.0, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)'))+
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

p8.pp50 <- ggplot(median.20, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) + 
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

p9.pp50 <- ggplot(median.50, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

png('results/plots/measures_pp50.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7.pp50,p8.pp50,p9.pp50)
dev.off()

# pp80 cutoff
performance.vector.0 <- c(summary0.pp80[1,],summary0.pp80[2,],summary0.pp80[3,],summary0.pp80[4,],summary0.pp80[5,],summary0.pp80[6,])
performance.vector.20 <- c(summary20.pp80[1,],summary20.pp80[2,],summary20.pp80[3,],summary0.pp80[4,],summary0.pp80[5,],summary0.pp80[6,])
performance.vector.50 <- c(summary50.pp80[1,],summary50.pp80[2,],summary50.pp80[3,],summary0.pp80[4,],summary0.pp80[5,],summary0.pp80[6,])
perc.1 <- rep(as.numeric(rownames(method.means)),6)*100
measure <- rep(c(rep('Precision',19),rep('Recall',19),rep('f-measure',19)),2)
method.1 <- c(rep('Lasso',57),rep('Semi',57))
group <- interaction(measure,method.1)
median.0 <- data.frame(percent=perc.1,performance=performance.vector.0,measure=measure,method=method.1,group=group)
median.20 <- data.frame(percent=perc.1,performance=performance.vector.20,measure=measure,method=method.1,group=group)
median.50 <- data.frame(percent=perc.1,performance=performance.vector.50,measure=measure,method=method.1,group=group)

p7.pp80 <- ggplot(median.0, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

p8.pp80 <- ggplot(median.20, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

p9.pp80 <- ggplot(median.50, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)')) +  
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))

png('results/plots/measures_pp80.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7.pp80,p8.pp80,p9.pp80)
dev.off()

# pp95 cutoff
performance.vector.0 <- c(summary0.pp95[1,],summary0.pp95[2,],summary0.pp95[3,],summary0.pp95[4,],summary0.pp95[5,],summary0.pp95[6,])
performance.vector.20 <- c(summary20.pp95[1,],summary20.pp95[2,],summary20.pp95[3,],summary0.pp95[4,],summary0.pp95[5,],summary0.pp95[6,])
performance.vector.50 <- c(summary50.pp95[1,],summary50.pp95[2,],summary50.pp95[3,],summary0.pp95[4,],summary0.pp95[5,],summary0.pp95[6,])
perc.1 <- rep(as.numeric(rownames(method.means)),6)*100
measure <- rep(c(rep('Precision',19),rep('Recall',19),rep('f-measure',19)),2)
method.1 <- c(rep('Lasso',57),rep('Semi',57))
group <- interaction(measure,method.1)
median.0 <- data.frame(percent=perc.1,performance=performance.vector.0,Measure=measure,Method=method.1,group=group)
median.20 <- data.frame(percent=perc.1,performance=performance.vector.20,Measure=measure,Method=method.1,group=group)
median.50 <- data.frame(percent=perc.1,performance=performance.vector.50,Measure=measure,Method=method.1,group=group)

p7.pp95 <- ggplot(median.0, aes(x=percent,y=performance,group=group,colour=Method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)',linetype='Measure: ',colour='Method: ')) +  
  ylim(0,1) + xlim(1,5) + theme(legend.position = c(0.5,0.8), legend.direction = "horizontal",legend.key.size = unit(1.5, "cm"),legend.text = element_text(size = 18),legend.title = element_text(size = 20),legend.key=element_blank()) + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=16,face = 'bold'), axis.title=element_text(size=16,face = 'bold'))

p8.pp95 <- ggplot(median.20, aes(x=percent,y=performance,group=group,colour=Method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)',linetype='Measure: ',colour='Method: ')) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = c(0.5,0.8), legend.direction = "horizontal",legend.key.size = unit(1.5, "cm"),legend.text = element_text(size = 18),legend.title = element_text(size = 20),legend.key=element_blank()) + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=16,face = 'bold'), axis.title=element_text(size=16,face = 'bold'))

p9.pp95 <- ggplot(median.50, aes(x=percent,y=performance,group=group,colour=Method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(y='Performance',x='Training Set (%)',linetype='Measure: ',colour='Method: ')) +  
  ylim(0,1) + xlim(1,5) + theme(legend.position = c(0.5,0.8), legend.direction = "horizontal",legend.key.size = unit(1.5, "cm"),legend.text = element_text(size = 18),legend.title = element_text(size = 20),legend.key=element_blank()) + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=16,face = 'bold'), axis.title=element_text(size=16,face = 'bold')) + labs(linetype='Measure',colour='Method')

png('results/plots/measures_pp95.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7.pp95,p8.pp95,p9.pp95)
dev.off()

# Save individual plots for precision, recall, and F-measure
png('results/plots/measures_median_0.png',height = 8, width = 8, res = 300, units = 'in')
p7
dev.off()

png('results/plots/measures_pp50_0.png',height = 8, width = 8, res = 300, units = 'in')
p7.pp50
dev.off()

png('results/plots/measures_pp80_0.png',height = 8, width = 8, res = 300, units = 'in')
p7.pp80
dev.off()

png('results/plots/measures_pp95_0.png',height = 8, width = 8, res = 300, units = 'in')
p7.pp95
dev.off()

png('results/plots/measures_median_20.png',height = 8, width = 8, res = 300, units = 'in')
p8
dev.off()

png('results/plots/measures_pp50_20.png',height = 8, width = 8, res = 300, units = 'in')
p8.pp50
dev.off()

png('results/plots/measures_pp80_20.png',height = 8, width = 8, res = 300, units = 'in')
p8.pp80
dev.off()

png('results/plots/measures_pp95_20.png',height = 8, width = 8, res = 300, units = 'in')
p8.pp95
dev.off()

png('results/plots/measures_median_50.png',height = 8, width = 8, res = 300, units = 'in')
p9
dev.off()

png('results/plots/measures_pp50_50.png',height = 8, width = 8, res = 300, units = 'in')
p9.pp50
dev.off()

png('results/plots/measures_pp80_50.png',height = 8, width = 8, res = 300, units = 'in')
p9.pp80
dev.off()

png('results/plots/measures_pp95_50.png',height = 8, width = 8, res = 300, units = 'in')
p9.pp95
dev.off()

# Midrange cutoff
performance.vector.0 <- c(summary0.midrange[1,],summary0.midrange[2,],summary0.midrange[3,],summary0.midrange[4,],summary0.midrange[5,],summary0.midrange[6,])
performance.vector.20 <- c(summary20.midrange[1,],summary20.midrange[2,],summary20.midrange[3,],summary0.midrange[4,],summary0.midrange[5,],summary0.midrange[6,])
performance.vector.50 <- c(summary50.midrange[1,],summary50.midrange[2,],summary50.midrange[3,],summary0.midrange[4,],summary0.midrange[5,],summary0.midrange[6,])
perc.1 <- rep(as.numeric(rownames(method.means)),6)*100
measure <- rep(c(rep('Precision',19),rep('Recall',19),rep('f-measure',19)),2)
method.1 <- c(rep('Lasso',57),rep('Semi',57))
group <- interaction(measure,method.1)
median.0 <- data.frame(percent=perc.1,performance=performance.vector.0,measure=measure,method=method.1,group=group)
median.20 <- data.frame(percent=perc.1,performance=performance.vector.20,measure=measure,method=method.1,group=group)
median.50 <- data.frame(percent=perc.1,performance=performance.vector.50,measure=measure,method=method.1,group=group)

p10 <- ggplot(median.0, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(title='Midrange Cutoff',y='Performance',x='Training Set (%)')) +  theme(plot.title=element_text(family="Times", face="bold", size=10)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed'))
p11 <- ggplot(median.20, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(title='Midrange Cutoff',y='Performance',x='Training Set (%)')) +  theme(plot.title=element_text(family="Times", face="bold", size=10)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed'))
p12 <- ggplot(median.50, aes(x=percent,y=performance,group=group,colour=method)) + geom_line(size=1,aes(linetype=measure)) + 
  labs(list(title='Midrange Cutoff',y='Performance',x='Training Set (%)')) +  theme(plot.title=element_text(family="Times", face="bold", size=10)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed'))

png('results/plots/measures_midrange.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p10,p11,p12)
dev.off()

# Transpose median, pp50, pp80, and pp95 cutoff

png('results/plots/cutoff0.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7,p7.pp50,p7.pp80,p7.pp95,ncol=1)
dev.off()

png('results/plots/cutoff20.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p8,p8.pp50,p8.pp80,p8.pp95,ncol=1)
dev.off()

png('results/plots/cutoff50.png',height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p9,p9.pp50,p9.pp80,p9.pp95,ncol=1)
dev.off()


# Table for number of positive predicted values
# LASSO
num.pos.pred.sup.0 <- round(cbind(apply(sup.pos.med.0,1,mean),
  apply(sup.pos.midrange.0,1,mean),
  apply(sup.pos.pp50.0,1,mean),
  apply(sup.pos.pp80.0,1,mean),
  apply(sup.pos.pp95.0,1,mean)),0)

colnames(num.pos.pred.sup.0) <- c('Median','Midrange','PP50','PP80','PP95')
rownames(num.pos.pred.sup.0) <- p*100

write.csv(num.pos.pred.sup.0,"results/tables/num_pos_pred_sup_0.csv")

num.pos.pred.sup.20 <- round(cbind(apply(sup.pos.med.20,1,mean),
                                  apply(sup.pos.midrange.20,1,mean),
                                  apply(sup.pos.pp50.20,1,mean),
                                  apply(sup.pos.pp80.20,1,mean),
                                  apply(sup.pos.pp95.20,1,mean)),0)

colnames(num.pos.pred.sup.20) <- c('Median','Midrange','PP50','PP80','PP95')
rownames(num.pos.pred.sup.20) <- p*100

write.csv(num.pos.pred.sup.20,"results/tables/num_pos_pred_sup_20.csv")

num.pos.pred.sup.50 <- round(cbind(apply(sup.pos.med.50,1,mean),
                                  apply(sup.pos.midrange.50,1,mean),
                                  apply(sup.pos.pp50.50,1,mean),
                                  apply(sup.pos.pp80.50,1,mean),
                                  apply(sup.pos.pp95.50,1,mean)),0)

colnames(num.pos.pred.sup.50) <- c('Median','Midrange','PP50','PP80','PP95')
rownames(num.pos.pred.sup.50) <- p*100

write.csv(num.pos.pred.sup.50,"results/tables/num_pos_pred_sup_50.csv")

# Semi-supervised
num.pos.pred.semi <- round(cbind(apply(semi.pos.med,1,mean),
  apply(semi.pos.pp50,1,mean),
  apply(semi.pos.pp80,1,mean),
  apply(semi.pos.pp95,1,mean)),0)
colnames(num.pos.pred.semi) <- c('Median','Midrange','PP50','PP80','PP95')
rownames(num.pos.pred.semi) <- p*100

write.csv(num.pos.pred.semi,"results/tables/num_pos_pred_semi.csv")

# save.image("saved_images/6D_boxplots_unbalanced")

