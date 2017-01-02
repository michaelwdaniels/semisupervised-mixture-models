# ---------------------------------------------------------------------------- #
# Script 0
# Authors: Rani Powers and Daniel Dvorkin
#
# Downloads the Seringhaus 2006 data from the Gerstein lab website, annotates
# each feature with appropriate group label (i.e. 'binary' or 'integer'), and 
# fits marginal models for each feature. Saves cleaned feature table and 
# marginal model fit report in results/tables directory.
# ---------------------------------------------------------------------------- #

# Get constants and helper functions
source("R/includes.R")
cat("SCRIPT 0: Download S. cerevisiae data, annotate features and fit marginals\n")

# --------------------------- GET CEREVISIAE DATA ---------------------------- #
# Download the Seringhaus 2006 data into the data/seringhaus/ directory
seringhaus.arff = "data/seringhaus/cerevisiae_ALL_noNaN.arff"
seringhaus.csv = "data/seringhaus/cerevisiae_ALL_noNaN.csv"
  
if (!file.exists(seringhaus.arff)){
  cat("\tDownloading Seringhaus S. cerevisiae data\n")
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
cat("\tFormatting S. cerevisiae data\n")
arffData = read.arff(seringhaus.arff)
featureNames = names(arffData)
featureNames = tolower(gsub("-", "_", featureNames))

# Read numeric feature values and break out into groups
features = read.table(seringhaus.csv, sep=",", as.is = TRUE)
names(features) = c("id", featureNames)

# Save formatted feature table
write.table(features, "data/seringhaus/cleaned_cerevisiae_features.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Extract gene IDs and essential gene identification
geneId = features$id
features$id = NULL
isEssential = as.logical(features$sgd_ess)
features$sgd_ess = NULL

# Remove non-sequence-dervived features, annotate, and adjust data that might 
# cause problems with log or logit transformation
features = mlapply(ALL_FEATURES, function(sdf)
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
names(features) = ALL_FEATURES

# --------------------------- FIT MARGINAL MODELS ---------------------------- #
# Marginal model fitting
cat("\tFitting marginal models\n")
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
#print(msrept)
msrept$feature = row.names(msrept)

# Save report on marginal models
write.table(msrept[,c(ncol(msrept), 1:(ncol(msrept)-1))], 
            "results/tables/cerevisiae_marginal_model_performance.txt",
            sep = '\t', row.names = F, quote = F)

cat("\tMarginal model stats saved in the 'results/tables/' directory\n")

# ------------------------- GET ESSENTIALITY DATA ---------------------------- #
# Read Saccharomyces Genome Database data for gene names and year discovered
sgd.ess = read.table("data/SGD/yeast_cerevisiae_essential_genes.tsv",
                     as.is=T, header=F, sep="\t")
sgd.ess = subset(sgd.ess, V2=="ORF") # open reading frames only
refs = sgd.ess$V5
pmid = as.integer(str_extract(refs, "[0-9]+"))
sgdref = str_extract(refs, "S[0-9]+")

# Keep only the desired columns
sgd.ess = sgd.ess[,c(1, 3, 9, 6)]
names(sgd.ess) = c("ensembl", "gene", "strain", "evidence")
sgd.ess$pmid = pmid
sgd.ess$sgdref = sgdref

# Get publications
sgd.pub = read.csv("Data/SGD/yeast_genome_database_publications.csv",
                   as.is=T, header=T)
names(sgd.pub) = c("author1", "pmid", "sgdref", "year")

# Merge essentiality and publication data, individual tables for two strains
sgd.big = subset(merge(sgd.ess, sgd.pub, all.x=TRUE), !is.na(year))
sgd288 = subset(sgd.big, strain=="S288C")
sgd1278 = subset(sgd.big, strain=="Sigma1278b")

# Save relevant info for analysis
preproc = list(features=features, essential=isEssential, summary=msrept,
               ensembl=geneId, sgd.big=sgd.big, sgd288=sgd288, sgd1278=sgd1278)
save(preproc, file="data/Rdata/cerevisiae_preprocessed.RData")
cat("\tPreprocessed S. cerevisiae features saved in the 'data/Rdata/ directory\n")