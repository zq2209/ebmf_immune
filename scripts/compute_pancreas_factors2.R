# This script should be run after running compute_pancreas_factors.R.
#
# sinteractive --mem=16G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(NNLM)
library(flashier)
library(fastTopics)
load("../data/pancreas.RData")
load("../output/pancreas_factors.RData")
set.seed(1)

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Remove genes with very low variance in expression.
x <- sparseMatrixStats::colSds(Y)
j <- which(x > 0.01)
Y <- Y[,j]

# Set alower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# (5) Fit an NMF with *cross-cutting factors* using flashier ("NMF-CC").
k <- 23
batch_topics <- c(2:5,7:8,20)
other_topics <- setdiff(1:k,batch_topics)
timings   <- list(fl_nmf_cc = 0)
t0        <- proc.time()
fl_nmf_cc <- flash_init(Y,var_type = 2,S = s1)
for (i in batch_topics) {
  fl_nmf_cc <- flash_factors_init(fl_nmf_cc,
                                  list(fl_nmf_ldf$L[,i,drop = FALSE],
                                       fl_nmf_ldf$F[,i,drop = FALSE]),
                                  c(ebnm_point_exponential,
                                    ebnm_point_normal))
}
for (i in other_topics) {
  fl_nmf_cc <- flash_factors_init(fl_nmf_cc,
                                  list(fl_nmf_ldf$L[,i,drop = FALSE],
                                       fl_nmf_ldf$F[,i,drop = FALSE]),
                                  ebnm_point_exponential)
}
fl_nmf_cc <- flash_backfit(fl_nmf_cc,extrapolate=FALSE,maxiter=100,verbose=3)
fl_nmf_cc <- flash_backfit(fl_nmf_cc,extrapolate=TRUE,maxiter=100,verbose=3)
t1 <- proc.time()
timings$fl_nmf_cc <- t1 - t0
print(timings$fl_nmf_cc)

# Save the model fit to an .Rdata file.
fl_nmf_cc_ldf <- ldf(fl_nmf_cc,type = "i")
session_info <- sessionInfo()
save(list = c("fl_nmf_cc_ldf","timings","session_info"),
     file = "pancreas_factors2.RData")
resaveRdaFiles("pancreas_factors2.RData")
