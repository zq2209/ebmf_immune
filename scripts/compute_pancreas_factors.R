# Analyze the pancreas data using various matrix factorization
# methods (flashier, fastTopics, etc)
#
# sinteractive --mem=16G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(NNLM)
library(fastTopics)
library(flashier)
load("../data/pancreas.RData")
set.seed(1)

# Remove genes that are expressed in fewer than 10 cells.
x      <- colSums(counts > 0)
j      <- which(x > 9)
counts <- counts[,j]

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- MatrixExtra::mapSparse(counts/(a*s),log1p)

# Remove genes with very low variance in expression (based on Y).
x <- sparseMatrixStats::colSds(Y)
j <- which(x > 0.01)
Y <- Y[,j]
counts <- counts[,j]

# Set a lower bound on the variances.
n  <- nrow(counts)
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Set up the "timings" data structure.
timings <- list(nmf        = 0,
                fl_nmf     = 0,
                fl_snmf    = 0,
                fasttopics = 0)

# (1) Fit an NMF using NNLM.
# I use k = 23 to match the flash() call immediately below.
Y_dense <- as.matrix(Y)
t0  <- proc.time()
nmf <- nnmf(Y_dense,k = 23,loss = "mse",method = "scd",
            max.iter = 200,verbose = 2,n.threads = 8)
t1  <- proc.time()
timings$nmf <- t1 - t0
print(timings$nmf)

# (2) Fit an NMF using flashier.
# I initially set var_type = 0 to increase the number of
# factors discovered.
t0 <- proc.time()
fl0 <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 0,
             greedy_Kmax = 40,nullcheck = FALSE,backfit = FALSE,
             verbose = 3)
fl_nmf <- flash_init(Y,var_type = 2,S = s1)
fl_nmf <- flash_factors_init(fl_nmf,fl0,ebnm_point_exponential)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 3)
t1 <- proc.time()
timings$fl_nmf <- t1 - t0
print(timings$fl_nmf)

# (3) Fit a semi-NMF using flashier.
# Here I set greedy_Kmax = 23 to align with the NMF flashier fit
# immediately above.
t0 <- proc.time()
fl0 <- flash(Y,ebnm_fn = c(ebnm_point_exponential,ebnm_point_normal),
             var_type = 0,greedy_Kmax = 23,nullcheck = FALSE,
             backfit = FALSE,verbose = 3)
fl_snmf <- flash_init(Y,var_type = 2,S = s1)
fl_snmf <- flash_factors_init(fl_snmf,fl0,
                              c(ebnm_point_exponential,
                                ebnm_point_normal))
fl_snmf <- flash_backfit(fl_snmf,extrapolate = FALSE,maxiter = 100,verbose = 3)
fl_snmf <- flash_backfit(fl_snmf,extrapolate = TRUE,maxiter = 100,verbose = 3)
t1 <- proc.time()
timings$fl_snmf <- t1 - t0
print(timings$fl_snmf)

# (4) Fit a Poisson NMF using fastTopics.
# I set k = 23 to align with the NMF flashier fit.
t0 <- proc.time()
pnmf0 <- fit_poisson_nmf(counts,k = 23,numiter = 100,method = "em",
                        control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                        init.method = "random",verbose = "detailed")
pnmf <- fit_poisson_nmf(counts,fit0 = pnmf0,numiter = 100,method = "scd",
                        control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                        verbose = "detailed")
t1 <- proc.time()
timings$fasttopics <- t1 - t0
print(timings$fasttopics)

# Save the model fits to an .Rdata file.
fl_nmf_ldf   <- ldf(fl_nmf,type = "i")
fl_snmf_ldf  <- ldf(fl_snmf,type = "i")
session_info <- sessionInfo()
save(list = c("nmf","fl_nmf_ldf","fl_snmf_ldf","pnmf","timings",
              "session_info"),
     file = "pancreas_factors.RData")
resaveRdaFiles("pancreas_factors.RData")
