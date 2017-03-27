### VISUALISATION ###
### This script provides three functions to generate principle component analysis (PCA) (from vcf-files), allele frequency spectra (AFS) (from vcftools outputs) and pairwise comparsion (see bash script 'pairwise_comparisons.sh' for details) plots, respectively.

library(data.table)
library(gdsfmt)
library(SNPRelate)

### PCA ###
# the customPCA() function requires packages 'gdsfmt' and 'SNPRelate' and expects merged vcf-files (usually one line of code using the PERL script 'vcf-merge' that comes with vcftools)

customPCA <- function(file, popFile) {
  set.seed(10)
  
  gds <- paste(file, '.gds', sep = '')
  snpgdsVCF2GDS(file, gds, method = "biallelic.only")
  #snpgdsSummary(gds)

  genofile <- openfn.gds(gds)

  snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2,
     autosome.only = FALSE)

  snpset.id <- unlist(snpset)

  pca <- snpgdsPCA(genofile, snp.id = snpset.id,
     autosome.only = FALSE)

  pc.percent <- pca$varprop*100

  sample.id <- read.gdsn(index.gdsn(genofile,
     "sample.id"))

  populations <- scan(popFile, what=character())

  tab <- data.frame(sample.id = pca$sample.id,
     pop = factor(populations)[match(pca$sample.id,
     sample.id)], EV1 = pca$eigenvect[,1],
     EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

  plot(tab$EV2, tab$EV1, col = as.integer(tab$pop), xlab = "eigenvector 2", ylab = "eigenvector 1", pch = 16)
  text(tab$EV2, tab$EV1, labels = tab$pop)  
  legend('bottomleft', legend = levels(tab$pop), pch = 16, col = 1:nlevels(tab$pop), bty = 'n')

  # output
  dev.copy(pdf, paste(file, '.pdf', sep = ''))
  dev.off()

  # done
  print('Output of PCA plot. Done.')
}


### ALLELE FREQUENCY SPECTRUM ###
system('vcftools --freq --out overall --gzvcf file.vcf.gz')
# the customAFS() function requires input from vcftools --freq calculations and package 'data.table', converts input into numeric values and plots a simple histogram

customAFS <- function(file) {
  freq <- fread(file)
 
  # convert to numeric
  newFreq <- sapply(freq$V6, function(x) as.numeric(strsplit(x, ':')[[1]][2]))
  
  # plot
  hist(newFreq, col = 4, xaxs = 'i', yaxs = 'i', xlab = 'Allele frequency', ylab = 'Counts', main = 'Allele frequency spectrum')

  # ouptut
  dev.copy(pdf, paste(file, '.pdf', sep = ''))
  dev.off()

  # done
  print('Output of AFS plot. Done.')
}


### HISTOGRAMS OF PAIRWISE DIFFERENCES ###
system('bash -x pairwise_differences.sh')
PD <- fread('PD.txt')
perPD <- PD$V4






