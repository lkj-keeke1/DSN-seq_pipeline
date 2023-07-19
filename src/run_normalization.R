#!/usr/bin/env Rscript


# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("Calculate FPKM/TPM/TMM based on featureCounts output(linux)")

# Add command line arguments
p <- add_argument(p, "featurecounts", help="input: featureCounts output table", type="character")
p <- add_argument(p, "outpre", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

fc_file <- argv$featurecounts
outFilePref <- argv$outpre

# libraries ---------------------------------------------------------------
library(tidyverse)
library(edgeR)


out_count_matrix_path <- paste(outFilePref, '.count_matrix', sep = '')
out_fpkm_matrix_path <- paste(outFilePref, '.fpkm_matrix', sep = '')
out_tpm_matrix_path <- paste(outFilePref, '.tpm_matrix', sep = '')
out_tmm_matrix_path <- paste(outFilePref, '.tmm_matrix', sep = '')
out_tmm_fpkm_matrix_path <- paste(outFilePref, '.tmm_fpkm_matrix', sep = '')
out_tmm_info_path <- paste(outFilePref, '.tmm_info', sep = '')

# Read featureCounts output matrix
fc_list <- read.table(file = fc_file, header = T, comment.char = "#") %>% 
  column_to_rownames(var = "Geneid")

# Split expression matrix and annotation information
 # Expression matrix
gene_exp <- fc_list[,6:ncol(fc_list)]
 # Gene annotation 
gene_ano <- fc_list[,1:5]


# fpkm calculation
 # Method 1: direct calculation
#nor_len <- gene_exp / gene_ano$Length
#fpkm <- t(t(nor_len)/colSums(gene_exp)) * 1e9

 # Method 2: edgeR package
dge_list <-  DGEList(counts=gene_exp, genes=gene_ano, group = colnames(gene_exp))
fpkm <-  rpkm(dge_list)


# tpm calculation 
 # Method 1: direct calculation
#nor_len <- gene_exp / gene_ano$Length
#tpm <- t(t(nor_len)/colSums(nor_len)) * 1e6

 # Method 2: from fpkm
tpm <- t(t(fpkm)/colSums(fpkm)) * 1e6


# TMM normalization (between samples normalization)
dge_list <- calcNormFactors(dge_list)
tmm <- cpm(dge_list)
tmm_fpkm <-  rpkm(dge_list)




# output

 # Expression matrixes
gene_exp <- rownames_to_column(as.data.frame(gene_exp), var = "Geneid")
fpkm <- rownames_to_column(as.data.frame(fpkm), var = "Geneid")
tpm <- rownames_to_column(as.data.frame(tpm), var = "Geneid")
tmm <- rownames_to_column(as.data.frame(tmm), var = "Geneid")
tmm_fpkm <- rownames_to_column(as.data.frame(tmm_fpkm), var = "Geneid")

write.table(x = gene_exp,
            file = out_count_matrix_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(x = fpkm,
            file = out_fpkm_matrix_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(x = tpm,
            file = out_tpm_matrix_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(x = tmm,
            file = out_tmm_matrix_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(x = tmm_fpkm,
            file = out_tmm_fpkm_matrix_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

 # TMM normalization factor information
dge_list$samples$eff.lib.size = dge_list$samples$lib.size * dge_list$samples$norm.factors
write.table(dge_list$samples,
            file = out_tmm_info_path,
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

