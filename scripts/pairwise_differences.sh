#!/bin/bash
### CALCULATES PAIRWISE DIFFERENCES BETWEEN VCF FILES.

### SORTING and FILTERING for biallelic SNPs in original SNPs - THIS NEEDS TO BE DONE FOR BOTH POPULATIONS ### 

# create list of all vcf-files
ls *vcf.gz > list_individuals

for ind in `cat list_individuals`; do
  # sorts vcf-file according to position
  vcf-sort "$ind" | bgzip -c > "$ind"_sorted.vcf.gz

  # only keep biallelic SNPs
  bcftools view -m2 -M2 -v snps "$ind"_sorted.vcf.gz | bgzip -c > "$ind"_sorted_biall.vcf.gz
done 

# clean-up
rm list_individuals
rm *sorted.vcf.gz


### PAIRWISE COMPARISON of original SNPs and SNP calls - this needs to be done separately for all three SNP callers, i.e. 'FreeB' with 'GATK' and 'Samtools' ###

# create a list of the sorted vcf-files
ls *biall* > list_PopA_biallelic

# for population A: ind1 corresponds the known real SNPs in a given individual (so the first loop goes over all 50 individuals), ind2 to SNP calls in that individual (so the second loop goes over all the SNP calls)
for ind1 in `cat popA_inds.txt`; do
  for ind2 in `cat list_PopA_biallelic`; do
    # collect ID of the individuals
    id1=$(echo "$ind1" | awk -F '[/d]' '{print $2}')
    id2=$(echo "$ind2" | awk -F '[/.]' '{print $1}' | awk -F '[/d]' '{print $2}')

    # only compare original SNPs and SNP calls in the same individuals
    if [ "$id1" == "$id2" ]
    then
      vcftools --gzvcf PopA_FreeB.vcf.gz --indv "$ind1" --gzdiff PopA_rSNP/"$ind2" --diff-site-discordance --out ${ind1}_${ind2}_diff
    fi
  done
done
# clean-up
rm list_PopA_biallelic

### COLLECT RESULTS from log files (check number of log files, must equal number of individuals!) ###
ls *log* > list_log

###  false positive rate
for log in `cat list_log`; do
  # grep number of different SNPs, i.e. SNPs that are different only present in the SNP calls but not in the simulated real data from log
  onlycalls=$(grep 'second file' "$log" | awk '{print $2}')
 
  # grep total number of SNPs from log
  total=$(grep 'possible' "$log" | awk '{print $4}')
  
  # generate ID for every individual for ouptut
  temp=$(grep '[--]out' "$log" | awk -F '[/.]' '{print $1, $6}'); inds=$(echo $temp | sed -e 's/gz_//g')
  
  # calculate false positive rate as proportion of different SNPs 
  PD=$(awk -v onlycalls="$onlycalls" -v total="$total" 'BEGIN {print onlycalls / total}'); echo $inds $PD >> PD.txt
done

# clean-up
rm list_log

ls *log* > list_log

### false negative rate - this is very similar to false positive rate but greps different numbers from the log
for log in `cat list_log`; do
  onlycalls=$(grep 'second file' "$log" | awk '{print $2}')

  # these are SNPs shared between the SNP calls and the simulated real data
  common=$(grep 'common' "$log" | awk '{print $2}')

  # this calculates the total number of variants from a vcf-file
  total=$(zgrep -v ^'#' "$file".log | wc -l)
  temp=$(grep '[--]out' "$log" | awk -F '[/.]' '{print $1, $6}'); inds=$(echo $temp | sed -e 's/gz_//g')
  
  # calculate false negative rate as total number of variants over number of false positives + shared SNPs, i.e. gives a rate of missing SNPs in the SNP calls
  PD=$(awk -v onlycalls="$onlycalls" -v common="$common" -v total="$total" 'BEGIN {print (total / (onlycalls + common))}'); echo $inds $PD >> PD_FN.txt
done

# clean-up
rm list_log


