# Fit semi-non-negative matrix factorizations (semi-NMFs) to the
# pancreas celseq2 data using flashier and gbcd.
#
# Note that I installed gbcd from the "form-YYT-option" branch:
# > remotes::install_github("stephenslab/gbcd@form-YYT-option")
# > packageVersion("gbcd")
# 0.2.4
library(tools)
library(Matrix)
library(ebnm)
library(flashier)
library(gbcd)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("../data/pancreas.RData")
set.seed(1)

# Select the CEL-seq2 data (Muraro et al, 2016).
# This should select 2,285 cells.
i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]

# Remove genes that are expressed in fewer than 10 cells.
x      <- colSums(counts > 0)
j      <- which(x > 9)
counts <- counts[,j]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# (1) Fit a semi-NMF of the *covariance matrix* using gbcd.
# Note that Kmax = 6 results in a factorization with at most 11 factors.
fl_cd <- fit_gbcd(Y,Kmax = 6,form_YYT = TRUE,
                  prior = flash_ebnm(prior_family = "point_exponential"),
                  maxiter1 = 100,maxiter2 = 100,maxiter3 = 100,
                  verbose = 3)

# (2) Fit a semi-NMF to Y using flashier, with at most 11 factors.
fl0 <- flash(Y,ebnm_fn = c(ebnm_point_exponential,ebnm_point_normal),
             var_type = 0,greedy_Kmax = 11,nullcheck = FALSE,
             backfit = FALSE,verbose = 3)
fl_snmf <- flash_init(Y,var_type = 2,S = s1)
fl_snmf <- flash_factors_init(fl_snmf,fl0,
                              c(ebnm_point_exponential,
                                ebnm_point_normal))
fl_snmf <- flash_backfit(fl_snmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_snmf <- flash_backfit(fl_snmf,extrapolate = TRUE,maxiter = 100,verbose = 3)

# Save the model fits to an .Rdata file.
fl_cd_ldf <- list(L = fl_cd$L,F = fl_cd$F$lfc)
fl_snmf_ldf  <- ldf(fl_snmf,type = "i")
session_info <- sessionInfo()
save(list = c("fl_cd_ldf","fl_snmf_ldf","session_info"),
     file = "pancreas_celseq2_snmf.RData")
resaveRdaFiles("pancreas_celseq2_snmf.RData")
