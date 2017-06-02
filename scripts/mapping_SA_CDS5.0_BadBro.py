#!/usr/bin/python

import os
import argparse
import glob
from os.path import basename
import re

####################
#   Informations   #
####################

courriel = "stephanie.arnoux@avignon.inra.fr"

parser=argparse.ArgumentParser(		#creation of a parser which gather the arguments of the program
                               description=""" program to map the reads of paired fastq files""",
                               epilog="""The input files must be in "fastq" format.\n
                                   contact: """+courriel+"""\n\n""")

parser.add_argument("-i",type=str,default="/data_nas1/wgs/solanaceae/Solution/raw_files/tomato/PAIRED_FILES/", 							help="input directory (default: %(default)s)")
parser.add_argument("-t",type=int,default=8,													help="threads (default: %(default)s)")
parser.add_argument("-mem",type=str,default="8g", 												help="maximum memory use (default: %(default)s)!")
parser.add_argument("-pr",type=str,default="LA", 												help="prefix (default: %(default)s)")
parser.add_argument("-ref",type=str,default="/data_nas1/wgs/solanaceae/Solution/ref_genomes/tomato/ITAG3.10_CDS.fasta", 		help="reference file (default: %(default)s)")
parser.add_argument("-refname",type=str,default="/data_nas1/wgs/solanaceae/Solution/ref_genomes/tomato/ITAG3.10_CDS",         help="reference file (default: %(default)s)")
parser.add_argument("-map",type=str,default="../mapping/tomatoCDS.3.1/", 								help="mapping directory (default: %(default)s)")
parser.add_argument("-rg",type=str,default="../mapping/tomatoCDS.3.1/rg_files/",                                 help="rg post mapping directory (default: %(default)s)")
parser.add_argument("-picard",type=str,	default="./picard-tools/picard.jar", 				help="picard location (default: %(default)s)")
parser.add_argument("-bwa",type=str,default="./bwa-0.7.12/bwa",help="bwa location(default: %(default)s)")
parser.add_argument("-samtools", 		type=str,		default="./samtools-0.1.19/samtools" , help="samtools location(default: %(default)s)")
parser.add_argument("-gatk",        type=str,       default="/usr/local/public/bi/GATK/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar" , help="gatk location(default: %(default)s)")
parser.add_argument("-vcftools",        type=str,       default="/usr/local/public/bi/vcftools/bin/vcftools" , help="vcftools location(default: %(default)s)")
parser.add_argument("-optionsHF",        type=str,       default='"QD<2.0||FS>60.0||MQ<40.0||ReadPosRankSum<-8.0"' , help="hard filtering option(default: %(default)s)")
parser.add_argument("-nameHF",        type=str,       default='"HardFilter"',  help="hard filtering name(default: %(default)s)")
parser.add_argument("-o",type=str,default="../vcf/tomatoCDS.3.1/",                    help="output directory (default: %(default)s)")

args=parser.parse_args()		#parsed arguments are gathered in args

################
#   Function   #
################
  
#Function try except  
def try_except(message,command):  
    try:  
        print(message+"\n")
        command = os.system(command)  
  
        if command == 0:  
            print(message+" done ...\n")  
        else:  
            print(message+" failed ...\n")  
    except OSError as e:  
        print("Execution failed ", e)             

#####################
#   Verifications   #
#####################

#Creating ../mapping/eggplant/statistics directory if it doesn"t exist
if not os.path.exists(args.map+"statistics"):
    try_except("the "+args.map+"statistics doesn\"t exists ... creating an empty directory","mkdir -p "+args.map+"statistics")

#Creating ../mapping/eggplant/rg_files directory if it doesn"t exist
if not os.path.exists(args.map+"rg_files"):
    try_except("the "+args.map+"rg_files doesn\"t exists ... creating an empty directory","mkdir -p "+args.map+"rg_files")

#Creating ../vcf/eggplant directory if it doesn"t exist
if not os.path.exists(args.o):
    try_except("the "+args.o+"  doesn\"t exists ... creating an empty directory","mkdir -p "+args.o)




#Creating index if one of the file of the index is missing
if os.path.isfile(args.ref):
    if not (os.path.isfile(args.ref+".amb") and os.path.isfile(args.ref+".ann") and os.path.isfile(args.ref+".bwt") and os.path.isfile(args.ref+".pac") and os.path.isfile(args.ref+".sa")):
    	try_except("creating bwa index files for "+args.ref+" ...",args.bwa+" index "+args.ref)
        if not os.path.isfile(args.ref+".fai") :
                try_except("creating samtools index files for "+args.ref+" ...","samtools faidx "+args.ref)
    elif not (os.path.isfile(args.refname+".dict")):
        try_except("generating picard index for "+args.refname+" ...","java_1.8 -Xmx"+args.mem+" -jar "+args.picard+" CreateSequenceDictionary R="+args.ref+"  O="+args.refname+".dict")

###################
#   Instructions  #
###################

#alignment
#try_except("paired end mapping for ...", "parallel --gnu -j"+str(args.t)+" '"+args.bwa+" mem -M -t8 "+args.ref+" {} $(echo {} | sed 's/1.fastq/2.fastq/') > "+args.map+"$(echo $(basename {}) | sed 's/1.fastq/.sam/')' ::: "+args.i+args.pr+"*1.fastq")

#conversion from sam to bam
#try_except("converting sam files to bam files ...","parallel --gnu -j"+str(args.t)+" 'samtools view -bS {} > $(echo {} | sed 's/.sam/.bam/')' ::: "+args.map+args.pr+"*.sam")

#sam removal
#try_except("removing sam files ...","parallel --gnu -j"+str(args.t)+" 'rm -f {}' ::: "+args.map+args.pr+"*.sam")

#sort bam
#try_except("sorting bam files ...","parallel --gnu -j"+str(args.t)+" '"+args.samtools+" sort {} "+args.map+"sorted_$(basename {} | sed 's/.bam//')' ::: "+args.map+args.pr+"*.bam")

#read groups replacement
#try_except("replacing reads groups for sorted bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.picard+" AddOrReplaceReadGroups I={} O=$(echo {} | sed 's/sorted_/rg_/') RGSM=$(echo $(basename {})) RGLB=Solution RGPL=illumina RGPU=none VALIDATION_STRINGENCY=LENIENT' ::: "+args.map+"sorted_*.bam")

#sorted bam removal
#try_except("removing sorted bam files...","parallel --gnu -j"+str(args.t)+" 'rm -f {}' ::: "+args.map+"sorted_*.bam")

#Move the rg_ files to the /rg_files/ folder
#try_except("moving rg bam files...","parallel --gnu -j"+str(args.t)+" 'mv -f {} "+args.map+"rg_files/$(echo $(basename {}))' ::: "+args.map+"rg_*.bam")

#bam index generation
#try_except("generating index for rg bam files ...","parallel --gnu -j"+str(args.t)+" '"+args.samtools+" index {}' ::: "+args.rg+"rg_*.bam")

#Mark duplicates
#try_except("generating metrics statistics for rg bam files and rmdup.bam files...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.picard+" MarkDuplicates I={} O=$(echo {} | sed 's/.bam/.rmdup.bam/') METRICS_FILE=$(echo {} | sed 's/.bam/.metrics/') MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ' ::: "+args.rg+"rg_*.bam")

#rmdup.bam index generation
#try_except("generating index for rmdup files...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.picard+" BuildBamIndex I={} ' ::: "+args.rg+"rg_"+args.pr+"*.rmdup.bam")

#idx statistics generation
#try_except("generating idx statistics for rg bam files ...","parallel --gnu -j"+str(args.t)+" 'samtools idxstats {} > "+args.map+"statistics/$(echo $(basename {}) | sed 's/rg_/idx_/' | sed 's/.bam/.txt/')' ::: "+args.rg+"rg_*rmdup.bam")

#flag statistics generation
#try_except("generating flag statistics for rg bam files ...","parallel --gnu -j"+str(args.t)+" 'samtools flagstat {} > "+args.map+"statistics/$(echo $(basename {}) | sed 's/rg_/flag_/' | sed 's/.bam/.txt/')' ::: "+args.rg+"rg_*rmdup.bam")

#intervals generation
try_except("generating intervals file for rg bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T RealignerTargetCreator -R "+args.ref+" -o $(echo {} | sed 's/.bam/.intervals/') -I {}' ::: "+args.rg+"rg_"+args.pr+"*.rmdup.bam")

#realigning rg bam files
try_except("realigning rg bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T IndelRealigner -R "+args.ref+" -targetIntervals $(echo {} | sed 's/.bam/.intervals/') -o "+os.getcwd()+"/$(echo {} | sed 's/.rmdup.bam/.rmdup.realn.bam/') -I {} -maxReads 100000' ::: "+args.rg+"rg_"+args.pr+"*.rmdup.bam")

#removing intervals
#try_except("removing intervals ...","parallel --gnu -j"+str(args.t)+" 'rm -f {}' ::: "+args.rg+"rg_*.intervals")

#removing rg bam files
#try_except("removing rg bam files ...","parallel --gnu -j"+str(args.t)+" 'rm -f {}' ::: "+args.rg+"rg_*.bam")

#removing rg bai files
#try_except("removing rg bai files ...","parallel --gnu -j"+str(args.t)+" 'rm -f {}' ::: "+args.rg+"rg_*.bai")

#calling variants for realigned bam files
try_except("calling variants for realigned bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T HaplotypeCaller -R "+args.ref+" -I {} -ERC GVCF   -variant_index_type LINEAR -variant_index_parameter 128000 -stand_emit_conf 20 -stand_call_conf 20 -o "+args.o+"$(echo $(basename {}) | sed 's/.realn.bam/.prerecal.gvcf/') ' ::: "+args.rg+"rg_"+args.pr+"*realn.bam")

#Generation of a list of files
try_except("Generation of a list of files","ls "+args.o+"*.prerecal.gvcf > "+args.o+args.pr+"_realn.list")

#Generation of a vcf really concervative for the BQSR
try_except("create vcf from gvcf files","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T GenotypeGVCFs -R "+args.ref+" -V {} -o "+args.o+"$(echo $(basename {}) | sed 's/.list/.Raw.vcf/') ' ::: "+args.o+args.pr+"_realn.list")

#hard filtering on the vcf file
try_except("create HARD filtered vcf","java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T VariantFiltration -R "+args.ref+" -V "+args.o+args.pr+"_realn.Raw.vcf --filterExpression "+args.optionsHF+" --filterName "+args.nameHF+" -o "+args.o+args.pr+"_realn.filtered.snps.vcf")

#MAF samtools filtering on the vcf file
try_except("create HARD+MAF filtered vcf",args.vcftools+" --vcf "+args.o+args.pr+"_realn.filtered.snps.vcf --out "+args.o+args.pr+"_realn.filtered.snps.vcf --remove-filtered-geno-all --maf 0.1 --recode ")  

#Base Recalibration
try_except("Base Recalibration ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T BaseRecalibrator -nct 12 -I {} --knownSites "+args.o+args.pr+"_realn.filtered.snps.vcf.recode.vcf -R "+args.ref+" -o $(echo {} | sed 's/.rmdup.realn.bam/.recal.data.table/')' ::: "+args.rg+"rg_*.rmdup.realn.bam")

#Print Reads (GATK Report Recalibration)
try_except("Generating printreads file for .realn.bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T PrintReads -nct 12 -I {} -BQSR $(echo {} | sed 's/.rmdup.realn.bam/.recal.data.table/') -R "+args.ref+" -o $(echo {} | sed 's/.rmdup.realn.bam/.recal.bam/')' ::: "+args.rg+"rg_*.rmdup.realn.bam")

#bam index generation
try_except("generating index for recal.bam files ...","parallel --gnu -j"+str(args.t)+" '"+args.samtools+" index {}' ::: "+args.rg+"*.recal.bam")

#calling variants for realigned bam files
try_except("calling variants for realigned bam files ...","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T HaplotypeCaller -R "+args.ref+" -I {} -ERC GVCF   -variant_index_type LINEAR -variant_index_parameter 128000 -stand_emit_conf 20 -stand_call_conf 20 -o "+args.o+"$(echo $(basename {}) | sed 's/.recal.bam/.3.1.gvcf/') ' ::: "+args.rg+"rg_"+args.pr+"*recal.bam")

#Generation of a list of final files :D
try_except("Generation of a list of files","ls "+args.o+"*.3.1.gvcf > "+args.o+args.pr+"_BaseCalibrated_ReadyToUse.list")

#Generation of gvcf of good quality
try_except("create vcf from gvcf files","parallel --gnu -j"+str(args.t)+" 'java_1.8 -Xmx"+args.mem+" -jar "+args.gatk+" -T GenotypeGVCFs -R "+args.ref+" -V {} -o "+args.o+"$(echo $(basename {}) | sed 's/_ReadyToUse.list/.vcf/') ' ::: "+args.o+args.pr+"_BaseCalibrated_ReadyToUse.list")

