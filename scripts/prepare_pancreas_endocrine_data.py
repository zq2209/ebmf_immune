# This script is adapted from the Jupyter notebook
# scRNA_seq_qc_preprocessing_clustering.ipynb found at
# https://github.com/theislab/pancreatic-endocrinogenesis
#
# To run this script, first download GSE132188_RAW.tar and extract the
# cellranger files from the GEO website, accession GSE132188, then put
# the files in the appropriate data subdirectories.
#
# Then run: python3 prepare_pancreas_endocrine_data.py
#
import numpy as np
import scipy as sci
import scanpy as sc
from anndata import AnnData
import anndata as ad

# Read the cellranger files for all four samples.
filename = "../data/GSE132188_RAW/E12_5_counts/mm10/matrix.mtx"
filename_genes = "../data/GSE132188_RAW/E12_5_counts/mm10/genes.tsv"
filename_barcodes = "../data/GSE132188_RAW/E12_5_counts/mm10/barcodes.tsv"

e125 = sc.read(filename).transpose()
e125.var_names = np.genfromtxt(filename_genes,dtype = str)[:, 1]
e125.obs_names = "e125-" + np.genfromtxt(filename_barcodes,dtype = str)

filename = "../data/GSE132188_RAW/E13_5_counts/mm10/matrix.mtx"
filename_genes = "../data/GSE132188_RAW/E13_5_counts/mm10/genes.tsv"
filename_barcodes = "../data/GSE132188_RAW/E13_5_counts/mm10/barcodes.tsv"

e135 = sc.read(filename).transpose()
e135.var_names = np.genfromtxt(filename_genes,dtype = str)[:, 1]
e135.obs_names = "e135-" + np.genfromtxt(filename_barcodes,dtype = str)

filename = "../data/GSE132188_RAW/E14_5_counts/mm10/matrix.mtx"
filename_genes = "../data/GSE132188_RAW/E14_5_counts/mm10/genes.tsv"
filename_barcodes = "../data/GSE132188_RAW/E14_5_counts/mm10/barcodes.tsv"

e145 = sc.read(filename).transpose()
e145.var_names = np.genfromtxt(filename_genes,dtype = str)[:, 1]
e145.obs_names = "e145-" + np.genfromtxt(filename_barcodes,dtype = str)

filename = "../data/GSE132188_RAW/E15_5_counts/mm10/matrix.mtx"
filename_genes = "../data/GSE132188_RAW/E15_5_counts/mm10/genes.tsv"
filename_barcodes = "../data/GSE132188_RAW/E15_5_counts/mm10/barcodes.tsv"

e155 = sc.read(filename).transpose()
e155.var_names = np.genfromtxt(filename_genes,dtype = str)[:, 1]
e155.obs_names = "e155-" + np.genfromtxt(filename_barcodes,dtype = str)

# Add dev. timepoint label for each sample.
e125.obs["day"] = "12.5"
e135.obs["day"] = "13.5"
e145.obs["day"] = "14.5"
e155.obs["day"] = "15.5"

# Create Concatenated anndata object for all timepoints.
e125.obs_names_make_unique()
e135.obs_names_make_unique()
e145.obs_names_make_unique()
e155.obs_names_make_unique()
e125.var_names_make_unique()
e135.var_names_make_unique()
e145.var_names_make_unique()
e155.var_names_make_unique()
alldays = ad.concat([e125, e135, e145, e155])

# Deleting individual day arrays.
del e125
del e135
del e145
del e155

# Quality control: calculate QC covariates for all anndata objects.
print(alldays.obs["day"].value_counts())

alldays.obs["n_counts"]   = alldays.X.sum(1)
alldays.obs["log_counts"] = np.log(alldays.obs["n_counts"])
alldays.obs["n_genes"]    = (alldays.X > 0).sum(1)

# Mitochondrial gene fraction.
mt_gene_mask = [gene.startswith("mt-") for gene in alldays.var_names]
mt_gene_index = np.where(mt_gene_mask)[0]
alldays.obs["mt_frac"] = alldays.X[:,mt_gene_index].sum(1) / alldays.X.sum(1)

# Filter cells according to identified QC thresholds.
print("Total number of cells: {:d}".format(alldays.n_obs))
alldays = alldays[alldays.obs["mt_frac"] < 0.2]
print("Number of cells after MT filter: {:d}".format(alldays.n_obs))

sc.pp.filter_cells(alldays,min_genes = 1200)
print("Number of cells after gene filter: {:d}".format(alldays.n_obs))
print("Total number of genes: {:d}".format(alldays.n_vars))

sc.pp.filter_genes(alldays,min_cells = 20)
print("Number of genes after cell filter: {:d}".format(alldays.n_vars))

# Write the filtered data.
alldays.write("pancreas_endocrine_alldays.h5ad")
