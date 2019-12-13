library(tidyverse)
library(qvalue)
library(svDialogs)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(enrichR)
library(matrixTests)

num_sample <- 11
num_exp <- 3
num_input <- 2
num_group <- 2
num_pool <- 2

# Normalize
x <- load_normalize(proteins_raw, peptides_raw, num_input, num_group, num_exp, num_pool)
y.prot <- log2_normalize(x$proteins_load_normalized[[1]], num_sample)
z.prot <- pool_normalize(y.prot, num_sample)
y.pep <- log2_normalize(x$peptides_load_normalized[[1]], num_sample)
z.pep <- pool_normalize(y.pep, num_sample)

# Relative Occupancy
m.pep <- occupancy(z.pep, z.prot)

# Univariate Stats
n.pep <- univariate(m.pep, num_sample, num_exp, num_pool)
n.prot <- univariate(z.prot, num_sample, num_exp, num_pool)

# Volcano Plot
o.pep <- select(n.pep, GN, `B / A Log2 Fold Change`, `B / A Adjusted P-Value`)
volcano(o.pep, p_val_cutoff = 0.15)

# Topn
topn(o.pep, topn = 25)

# Proteins vs. Peptides
p.pep <- select(n.pep, Accession, `B / A Log2 Fold Change`, `B / A Adjusted P-Value`)
p.prot <- select(n.prot, Accession, `B / A Log2 Fold Change`, `B / A Adjusted P-Value`)
protein(p.prot, p.pep)

# Enrich
enrich(o.pep, 100, p_cutoff = 1, enrichr_p_cutoff = 1)

