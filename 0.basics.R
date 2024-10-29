

# Imports/environment -----------------------------------------------------


library(ggplot2); library(reshape2); library(ggpubr); library(gridExtra); library(vcfR)
require(RColorBrewer); require(scales); library(stringr); library(betareg)
library(jtools); library(betareg); library(phytools); library(matrixStats);
library(GenomicRanges); library(rtracklayer); library(readxl); library(tidyverse)
library(lmerTest)

options(stringsAsFactors = F)
setwd('~/EPICC/') # Change this folder depending on the location of the downloaded files

theme_mypub <- function(base_size = 14,
                        base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

# Sample lists ------------------------------------------------------------

epicc.df <- read.delim('ListDNAPass.EPICC.txt')
samples <- unique(epicc.df$Patient) # patient ID list
msiList <- c('C552','C518','C516','C548','C536','C562')
adenomaList <- unique(epicc.df$Patient[epicc.df$Tissue=='Adenoma'])

rna.df <- read.delim('RNA/ListRNAPass.EPICC.txt')
atac.df <- read.delim('ATAC/ListATAC.EPICC.txt')

ip.samples.df <- read.delim('FFPE_samples/immunopanel_sample_key.txt')
patientEsc.df <- read.delim('Immune_escape/EPICC_escape_perPatient.manual.txt')


