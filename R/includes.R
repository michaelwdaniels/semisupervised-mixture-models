# --------------------------------- PACKAGES --------------------------------- #

# All packages needed for the full analysis
library(lcmix)
library(foreign)     # for read.arff
library(stringr)     
library(parallel)    
library(plyr)
library(RWeka)
library(ROCR)
library(reshape2)

# ------------------------------- "CONSTANTS" -------------------------------- #

# Number of cores to be used in parallel processes
NCORES = max(1, detectCores()-2) 

# S. cerevisiae features used in this project
SEQUENCE_DERIVED_FEATURES = c("mitochondria", "cytoplasm", "er", "nucleus", 
                              "vacuole", "other", "tm_helix", "l_aa", "gravy", 
                              "nc", "cai", "gc", "rare_aa_ratio", "close_stop_ratio")
ALL_FEATURES = c(SEQUENCE_DERIVED_FEATURES, "intron", "intxn_partners", "chromosome",
                 "in_how_many_of_5_proks_blast", "in_how_many_of_6_close_yeast_blast",
                 "blast_hits_in_yeast", "dovexpr", "chr_position")

# Define a group for each feature. Scripts will annotate features accordingly
# prior to fitting marginal and joint models
BINARY_FEATURES = namedVector("mitochondria", "cytoplasm", "er", "nucleus", 
                              "vacuole", "other", "intron")
INTEGER_FEATURES = namedVector("tm_helix", "l_aa", "intxn_partners", "chromosome",
                               "in_how_many_of_5_proks_blast",
                               "in_how_many_of_6_close_yeast_blast",
                               "blast_hits_in_yeast")
OPEN_FEATURES = namedVector("gravy", "nc", "dovexpr")
CLOSED_FEATURES = namedVector("cai", "gc", "rare_aa_ratio", "close_stop_ratio", 
                              "chr_position")

# Discrete as binary or integer; continuous on the open interval of the entire 
# real line, the half-open interval, and the closed interval [0,1]. Note that 
# the distributions of all the [0,1] features appear to be sufficiently 
# "gamma-like" that we can analyze them on the half-open interval.
DISCRETE_FEATURES = c(BINARY_FEATURES, INTEGER_FEATURES)
CONTINUOUS_FEATURES = c(OPEN_FEATURES, CLOSED_FEATURES)

# Baseline fits for marginal models
BASELINE_FIT_FAMILY = list(
  binary = "bernoulli",
  integer = "poisson",
  open = "normal",
  closed = "gamma")


# ------------------------------- FUNCTIONS ---------------------------------- #

binaryhist <- function(x, colors=1:length(x), ...){
  # Plot density histogram for the elements of list "x" given that such elements 
  # are binary (i.e., all elements of each element of x are 0 or 1)
  
  I = length(x)
  width = 0.1
  
  hists = lapply(1:I, function(i){
    res = hist(x[[i]], breaks=c(0, width, 1-width, 1), plot=FALSE)
    res$density = res$density / sum(res$density)
    res$intensities = res$density
    res$counts = res$density
    res$breaks = res$breaks + (i-1)*width    # shift x-axis values
    return(res)
  })
  
  # hists = lapply(x, hist, breaks=breaks, freq=FALSE, plot=FALSE)
  ylim = c(0, max(sapply(hists, function(h) h$counts))) * 1.05
  # leave room above bars so we don"t have highest bar all the way to the top, which can be deceptive
  plot(hists[[1]], xlim=c(0, 0.9+(I*width)), ylim=ylim,
       axes=FALSE, freq=FALSE, col=colors[[1]], ...)
  if(I > 1) for(i in 2:I) plot(hists[[i]], col=colors[[i]], add=TRUE)
  box()
  
  axis(2)
  axis(1, at=c(I*width/2, 1 + (I-2)*width/2), label=c(0,1))
}


mlapply <- function(X, FUN, ...)
  # Multicore utility wrapper for mclapply()
  mclapply(X, FUN, ..., mc.cores=NCORES)


msapply <- function (X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
  # Multicore version of sapply() depending on mlapply()
  # Code is taken almost verbatim from sapply()
  
  FUN <- match.fun(FUN)
  answer <- mlapply(X=X, FUN=FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

  
name_by_rownames <- function(x, name="name"){
  # Given a matrix or data frame "x" with non-null row names, return a data 
  # frame with a leading name column equal to the row names.
  
  if(!isMDF(x) || is.null(rownames(x)))
    stop("'x' must be a matrix or data frame with non-NULL row names")
  
  addon = data.frame(name=rownames(x))
  names(addon) = name
  cbind(addon, as.data.frame(x))
}
