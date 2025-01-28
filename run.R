library(dplyr)
library(Matrix)
library(limma)
library(pheatmap)
library(utils)
library(BiocGenerics)
library(gtools)
library(stats)
library(data.table)
library(DWLS)
library(BisqueRNA)
library(MuSiC)
library(quadprog)
library(SingleCellExperiment)
library(SCDC)
library(MAST)
library(deconvSeq)
library(BayesPrism)

setwd("path/to/your/folder") # change to your folder

source("utils_method.R")
source("utils_benchmark.R")
rm(P)

# Read data and metadata
data = readRDS(list.files(path="baron",pattern="rds",full.names=TRUE))
full_phenoData=read.table(list.files(path="baron",pattern="phenoData",full.names=TRUE),header=TRUE)
#---r

# Simulate bulk option
number_cells = 100; # Number of cells to be used to make the pseudo-bulk mixtures (multiple of 100)
number_mixtures = 10; # Number of mixtures generated 
# transformation = "none"
# normalization = "none"

source("benchmark.R")