# Analysis the pancreas celseq2 data using various matrix
# factorization methods (flashier, fastTopics, etc)
library(tools)
library(Matrix)
library(NNLM)
library(flashier)
library(fastTopics)
load("../data/pancreas.RData")
set.seed(1)

# All the matrix factorizations will have this many topics or factors.
k <- 9

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

# (1) Fit an NMF using NNLM.
Y_dense <- as.matrix(Y)
nmf <- nnmf(Y_dense,k = k,loss = "mse",method = "scd",
            max.iter = 200,verbose = 2,n.threads = 8)

# (2) Fit an NMF using fastTopics.
pnmf0 <- fit_poisson_nmf(counts,k = k,numiter = 100,method = "em",
                        control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                        init.method = "random",verbose = "detailed")
pnmf <- fit_poisson_nmf(counts,fit0 = pnmf0,numiter = 100,method = "scd",
                        control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                        verbose = "detailed")

# (3) Fit an NMF using flashier.
# I initially set var_type = 0 to increase the number of
# factors discovered.
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = k,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl_nmf <- flash_init(Y,var_type = 2,S = s1)
fl_nmf <- flash_factors_init(fl_nmf,fl0,ebnm_point_exponential)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 3)

# (4) Fit a semi-NMF using flashier.
fl0 <- flash(Y,ebnm_fn = c(ebnm_point_exponential,ebnm_point_normal),
             var_type = 0,greedy_Kmax = k,nullcheck = FALSE,
             backfit = FALSE,verbose = 3)
fl_snmf <- flash_init(Y,var_type = 2,S = s1)
fl_snmf <- flash_factors_init(fl_snmf,fl0,
                              c(ebnm_point_exponential,
                                ebnm_point_normal))
fl_snmf <- flash_backfit(fl_snmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_snmf <- flash_backfit(fl_snmf,extrapolate = TRUE,maxiter = 100,verbose = 3)

# Save the model fits to an .Rdata file.
fl_nmf_ldf   <- ldf(fl_nmf,type = "i")
fl_snmf_ldf  <- ldf(fl_snmf,type = "i")
session_info <- sessionInfo()
save(list = c("nmf","fl_nmf_ldf","fl_snmf_ldf","pnmf",
              "session_info"),
     file = "pancreas_celseq2_factors.RData")
resaveRdaFiles("pancreas_celseq2_factors.RData")
