# Semisupervised mixture models
Semi-supervised mixture models with the `lcmix` R package for integrative genomics analysis

## Scripts 1A and 1B
Semi-supervised and unsupervised mixture models on either sequence-derived features only (1A) or all features (1B). Train models on pre-2002 essential genes (+ training examples only, n = 64) and test on post-2002 essential and non-essential genes (+/- examples). Each script plots an ROC curve comparing the semi-supervised mixture model to the unsupervised mixture model.

## Scripts 2A and 2B
Semi-supervised and unsupervised mixture models on either sequence-derived features only (2A) or all features (2B). Train models on essential genes (+ training examples only, varying number of training examples) and test on remaining essential and non-essential genes (+/- examples). Each script plots an ROC curve comparing semi-supervised mixture models and unsupervised mixture models at varying numbers of training examples.

## Scripts 3A and 3B
Semi-supervised and unsupervised mixture models on either sequence-derived features only (3A) or all features (3B). Train models on essential and non-essential genes (+/- training examples, varying number) and test on remaining essential and non-essential genes (+/- examples). Each script plots an ROC curve comparing semi-supervised mixture models and unsupervised mixture models at varying numbers of training examples.

## includes.R
Loads required packages and some `lcmix` helper functions.

## analysis_helpers.R
Includes wrapper functions for the modeling and plotting used in all the scripts above.

## Scripts 6D
Runs simulation for supervised comparisons (LASSO, SVM, RandomForest). Time warning: The program takes my MacBook Pro around 36 hours to run.
