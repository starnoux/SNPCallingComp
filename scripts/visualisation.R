### VISUALISATION ###
### This script provides three functions to generate principle component analysis (PCA) (from vcf-files), allele frequency spectra (AFS) (from vcftools outputs) and pairwise comparsion (see bash script 'pairwise_comparisons.sh' for details) plots, respectively.

library(data.table)
library(gdsfmt)
library(SNPRelate)

### PCA ###
# the customPCA() function requires packages 'gdsfmt' and 'SNPRelate' and expects merged vcf-files (usually one line of code using the PERL script 'vcf-merge' that comes with vcftools). the function is heavily built on the population genetics tutorial (credit to Susanne P. Pfeifer).

customPCA <- function(file, popFile) {
  #LD pruning below introduces stochasticity - setting a random seed gives reproducible results
  set.seed(10)
  
  # create gds-file name and convert biallelic sites only from vcf to gds format - then open.
  gds <- paste(file, '.gds', sep = '')
  snpgdsVCF2GDS(file, gds, method = "biallelic.only")
  genofile <- openfn.gds(gds)

  # LD pruning removes one variant of SNP pairs in high LD
  snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2,
     autosome.only = FALSE)
  snpset.id <- unlist(snpset)

  # calculate PCA from the snpset
  pca <- snpgdsPCA(genofile, snp.id = snpset.id,
     autosome.only = FALSE)
  pc.percent <- pca$varprop*100

  # assign individuals and prepare labels
  sample.id <- read.gdsn(index.gdsn(genofile,
     "sample.id"))
  populations <- scan(popFile, what=character())

  # convert pca-results to data.frame
  tab <- data.frame(sample.id = pca$sample.id,
     pop = factor(populations)[match(pca$sample.id,
     sample.id)], EV1 = pca$eigenvect[,1],
     EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

  # plot PCA results
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


### BOXPLOTS OF PAIRWISE DIFFERENCES ###
system('bash -x pairwise_differences.sh')
# to generate the boxplots, the output of the script 'pairwise_differences.sh' and the pacakge 'data.table'' are needed 

### false positives
PDFreeB <- fread('PD_FreeB.txt')
PDSamtools <- fread('PD_Samtools.txt')
PDGATK <- fread('PD_GATK.txt')

perPDFreeB <- PDFreeB$V4
perPDSamtools <- PDSamtools$V4
perPDGATK <- PDGATK$V4

# with outlier
boxplot(perPDGATK, perPDSamtools, perPDFreeB, xaxs = NULL, col = '#990000', ylab = 'Ratio of false positives to total numer of variants')
axis(1, at = 1:3, c('GATK', 'Samtools', 'FreeB'))

# output
dev.copy(pdf, 'falsePositives_with_outlier.pdf')
dev.off()

#without outlier
boxplot(perPDGATK, perPDSamtools, perPDFreeB, xaxs = NULL, col = '#990000', ylab = 'Ratio of false positives to total numer of variants', ylim = c(0, 0.0014))
axis(1, at = 1:3, c('GATK', 'Samtools', 'FreeB'))
dev.copy(pdf, 'falsePositives_wo_outlier.pdf')
dev.off()

### false negatives
PDFNFreeB <- fread('PD_FN_FreeB.txt')
PDFNSamtools <- fread('PD_FN_Samtools.txt')
PDFNGATK <- fread('PD_FN_GATK.txt')

perPDFNFreeB <- PDFNFreeB$V4
perPDFNSamtools <- PDFNSamtools$V4
perPDFNGATK <- PDFNGATK$V4

boxplot(perPDFNGATK, perPDFNSamtools, perPDFNFreeB, xaxs = NULL, col = '#990000', ylab = 'Ratio of false negatives to total numer of variants')
axis(1, at = 1:3, c('GATK', 'Samtools', 'FreeB'))
dev.copy(pdf, 'falseNegatives.pdf')
dev.off()






