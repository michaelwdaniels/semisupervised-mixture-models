# ---------------------------------------------------------------------------- #
# Script 0
# Authors: Rani Powers and Daniel Dvorkin
#
# Downloads the Seringhaus 2006 data from the Gerstein lab website, annotates
# each feature with appropriate group label (i.e. 'binary' or 'integer'), and 
# fits marginal models for each feature. Saves cleaned feature table and 
# marginal model fit report in results/tables directory.
# ---------------------------------------------------------------------------- #

cat("SCRIPT 0: Download yeast data, annotate features and fit marginals\n")

# Get constants and helper functions
source("R/includes.R")

# Download the Seringhaus 2006 data into the data/seringhaus/ directory
cat("\tDownloading Seringhaus data\n")
seringhaus.arff = "data/seringhaus/cerevisiae_ALL_noNaN.arff"
seringhaus.csv = "data/seringhaus/cerevisiae_ALL_noNaN.csv"
  
if (!file.exists(seringhaus.arff)){
  download.file(url = "http://www.gersteinlab.org/proj/predess/data/Scerevisiae/Compiled/cerevisiae_ALL_noNaN.arff", 
                destfile = seringhaus.arff, method = "curl")
}
if (!file.exists(seringhaus.csv)){
  download.file(url = "http://www.gersteinlab.org/proj/predess/data/Scerevisiae/Compiled/cerevisiae_ALL_noNaN.csv", 
                destfile = seringhaus.csv, method = "curl")
}

# Read cerevisiae feature data from .arff file
# The only reason we actually do this is to get feature names, since it's easier 
# to get the actual data from cerevisiae_compiled_features.csv, which includes 
# gene IDs. Note that if for some reason we did choose to use the numerical data 
# from the .arff file, read.arff() doesn't have an "as.is" option, so we'd need 
# to convert factors to the corresponsing numeric values
cat("\tFormatting Seringhaus data\n")
arffData = read.arff(seringhaus.arff)
featureNames = names(arffData)
featureNames = tolower(gsub("-", "_", featureNames))

# Read numeric feature values and break out into groups
features = read.table(seringhaus.csv, sep=",", as.is = TRUE)
names(features) = c("id", featureNames)

# Save formatted feature table
write.table(features, "data/seringhaus/cleaned_Seringhaus_features.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Extract gene IDs and essential gene identification
geneId = features$id
features$id = NULL
isEssential = as.logical(features$sgd_ess)
features$sgd_ess = NULL

# Remove non-sequence-dervived features, annotate, and adjust data that might 
# cause problems with log or logit transformation
features = features[,SEQUENCE_DERIVED_FEATURES]
features = mlapply(SEQUENCE_DERIVED_FEATURES, function(sdf)
{
  res = features[[sdf]]
  
  if(sdf %in% BINARY_FEATURES)
    attr(res, "group") = "binary"
  if(sdf %in% INTEGER_FEATURES)
    attr(res, "group") = "integer"
  if(sdf %in% OPEN_FEATURES)
    attr(res, "group") = "open"
  if(sdf %in% CLOSED_FEATURES) {
    res[res < LC_EPS] = LC_EPS
    res[res > LC_REPS] = LC_REPS
    attr(res, "group") = "closed"
  }
  
  attr(res, "continuous") = (sdf %in% CONTINUOUS_FEATURES)
  attr(res, "family") = BASELINE_FIT_FAMILY[[attr(res, "group")]]
  
  return(res)
})


# Marginal model fitting
cat("\tFitting marginal models\n\n")
decincfit <- function(x, K, ...){
  # function to determine decreasing vs. increasing by essential genes
  fits = list(mixmod(x, K, decreasing = TRUE, ...),
              mixmod(x, K, decreasing = FALSE, ...))
  aucs = sapply(fits, rocauc, labels = isEssential)
  fits[[which.max(aucs)]]
}

features = mlapply(features, function(x){
  # basic fitting
  family = attr(x, "family")
  fits = lapply(2:3, function(K) decincfit(x, K, family = family))
  iclbics = sapply(fits, iclbic, map = TRUE)
  imax = max(iclbics)
  Km1 = which.max(iclbics)
  bestfit = fits[[Km1]]
  K = 1 + Km1 # because fits[[1]] corresponds to K=2, etc.
  attr(x, "K") = K
  
  # test for alternative distributions
  if(family == "poisson") {
    altfit = decincfit(x, K, family="nbinom")
    ialt = iclbic(altfit, map=TRUE)
    if(ialt > imax) {
      attr(x, "family") = "nbinom"
      bestfit = altfit
    }
  } else if(family == "normal") {
    altfit = decincfit(x, K, family="pvii")
    ialt = iclbic(altfit, map=TRUE)
    if(ialt > imax) {
      attr(x, "family") = "pvii"
      bestfit = altfit
    }
  }
  
  # save useful attributes
  attr(x, "rocauc") = rocauc(bestfit, isEssential)
  attr(x, "decreasing") = bestfit$decreasing
  attr(x, "bic") = bic(bestfit)
  attr(x, "iclbic") = iclbic(bestfit, map=TRUE)
  
  # send it back
  return(x)
})

# Order features by predictive value
raucs = sapply(features, function(x) attr(x, "rocauc"))
features = features[order(raucs, decreasing=TRUE)]

# Report on model selection
attnames = setdiff(names(attributes(features[[1]])), "continuous")
names(attnames) = attnames
msrept = as.data.frame(lapply(attnames, function(an)
  sapply(features, function(x) attr(x, an))),
  stringsAsFactors = F)
print(msrept)
msrept$feature = row.names(msrept)

# Save report on marginal models
write.table(msrept[,c(ncol(msrept), 1:(ncol(msrept)-1))], 
            "results/tables/yeast_marginal_model_performance.txt",
            sep = '\t', row.names = F, quote = F)

cat("\n\tMarginal model fitting complete\n")
cat("\tMarginal model stats saved in the 'results/tables/' directory")

