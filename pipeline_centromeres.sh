#!/bin/bash -l
#SBATCH --job-name=centromeres_pipeline_Norin61
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=3-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=6G   # Memory per core, use --mem= for memory per node
#SBATCH --output="%x_%j.o"
#SBATCH --error="%x_%j.e"
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,killable.q





###################################################################################################################################


## demultiplex CHIP002 data

# module load GLib
# INDEX='/bulk/jpoland/centromeres/CHIP002/'
# /homes/jpoland/tools/deML/src/deML -i index_sample.txt -f NEX0008_R1.fastq.gz -r NEX0008_R2.fastq.gz -if1 NEX0008_I1_001.fastq.gz -if2 NEX0008_I2_001.fastq.gz -o demultiplexed





###################################################################################################################################


## CHINESE SPRING V2  ##
## BUILD BOWTIE
# FASTA="/bulk/jpoland/genome/wheat/ChineseSpring/index/CS_refseqv2.fasta"

# /homes/jpoland/bowtie2-2.3.4.3/bowtie2-build  ${FASTA}  ${INDEX}	>> ${OUTPUT}

## BUILD HISAT
# /homes/jpoland/tools/hisat2-2.1.0/hisat2-build ${FASTA}  ${INDEX}	>> ${OUTPUT}


## set directory variables


###################################################################################################################################

## Chinese Spring ##


## data from Guo et al. 2016 ##


# NAME="cs.guo"
# 
# INDEX="/bulk/jpoland/genome/wheat/ChineseSpring/index/CS_refseqv1"
# # INDEX2="/bulk/jpoland/genome/wheat/ChineseSpring/index/CS_refseqv2c"
# 
# 
# CHIP1="/bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_SRR1686799/SRR1686799_1.fastq"
# CHIP2="/bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_SRR1686799/SRR1686799_2.fastq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/chipseq/chinesespring-guo_R1.trim.fastq"
# CHIP2TRIM="/bulk/jpoland/centromeres/chipseq/chinesespring-guo_R2.trim.fastq"




###################################################################################################################################


## Chinese Spring data from Koo ##


# NAME="cs.koo"
# 
# INDEX="/bulk/jpoland/genome/wheat/ChineseSpring/index/CS_refseqv1"
# # INDEX2="/bulk/jpoland/genome/wheat/ChineseSpring/index/CS_refseqv2c"
# 
# 
# CHIP1="/bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_KSU-Koo/KSU_AA-ChipSeqCENH3_NoIndex_L001_R1_001.fastq.gz"
# CHIP2="/bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_KSU-Koo/KSU_AA-ChipSeqCENH3_NoIndex_L001_R2_001.fastq.gz"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/chipseq/chinesespring-koo_R1.fastq.trim.gz"
# CHIP2TRIM="/bulk/jpoland/centromeres/chipseq/chinesespring-koo_R2.fastq.trim.gz"
# 




###################################################################################################################################




## JAGGER ## data from Koo

## set directory variables
# SEQUENCE="/bulk/jpoland/sequence/"
# FASTA="/bulk/jpoland/genome/wheat/Jagger/index/180529_Jagger_pseudomolecule_v1.1.fasta"
# GENOME="/bulk/jpoland/genome/wheat/Jagger/index/180529_Jagger_pseudomolecule_v1.1"
# 
# NAME="jagger"
# 
# INDEX="/bulk/jpoland/genome/wheat/Jagger/index/180529_Jagger_pseudomolecule_v1.1.fasta"
# 
# CHIP1="/bulk/jpoland/genome/wheat/Jagger/raw_data/CHIP-CENH3/Jagger-CENH3_R1.fastq.gz"
# CHIP2="/bulk/jpoland/genome/wheat/Jagger/raw_data/CHIP-CENH3/Jagger-CENH3_R2.fastq.gz"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/chipseq/jagger_R1.fastq.gz"
# CHIP2TRIM="/bulk/jpoland/centromeres/chipseq/jagger_R2.fastq.gz"



###################################################################################################################################

## trim adapters from CHIP0002 library ##
# 
# CHIP1="/bulk/jpoland/globus/psomagen/1909UNHS-0017/CHIP002/CHIP002_R1.fastq.gz"
# CHIP2="/bulk/jpoland/globus/psomagen/1909UNHS-0017/CHIP002/CHIP002_R2.fastq.gz"
# 
# CHIP1TRIM="/bulk/jpoland/globus/psomagen/1909UNHS-0017/CHIP002/CHIP002_R1.trim.fastq.gz"
# CHIP2TRIM="/bulk/jpoland/globus/psomagen/1909UNHS-0017/CHIP002/CHIP002_R2.trim.fastq.gz"
# 


###################################################################################################################################

## ArinaLrFor ## 

## BUILD HISAT
# FASTA="/bulk/jpoland/genome/wheat/ArinaLrFor/index/180808_Arina_pseudomolecules_v3.fasta"
# INDEX="/bulk/jpoland/genome/wheat/ArinaLrFor/index/1180808_Arina_pseudomolecules_v3"
#  
# /homes/jpoland/tools/hisat2-2.1.0/hisat2-build -p 8 ${FASTA}  ${INDEX}
# 

# 
# NAME="arina"
# 
# INDEX="/bulk/jpoland/genome/wheat/ArinaLrFor/index/1180808_Arina_pseudomolecules_v3"
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/ArinaLrFor-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/ArinaLrFor-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/ArinaLrFor-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/ArinaLrFor-CENH3_R2.trim.fq"
# 


###################################################################################################################################

## Julius ##

# # FASTA="/bulk/jpoland/genome/wheat/Julius/index/170807_julius_MAGIC3_pseudomolecules_v1.fasta"
# # INDEX="/bulk/jpoland/genome/wheat/Julius/index/170807_julius_MAGIC3_pseudomolecules_v1"
#  
# ## BUILD HISAT
# # /homes/jpoland/tools/hisat2-2.1.0/hisat2-build -p 8 ${FASTA}  ${INDEX}
# 

# 
# NAME="julius"
#  
# INDEX="/bulk/jpoland/genome/wheat/Julius/index/170807_julius_MAGIC3_pseudomolecules_v1"
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/Julius-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/Julius-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/Julius-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/Julius-CENH3_R2.trim.fq"
# 



###################################################################################################################################

## Lancer ##

# NAME="lancer"
# 
# INDEX="/bulk/jpoland/genome/wheat/Lancer/index/181120_lancer_pseudomolecule_v1.0"
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/Lancer-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/Lancer-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/Lancer-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/Lancer-CENH3_R2.trim.fq"
# 



###################################################################################################################################

## CDC Landmark ## 

# NAME="landmark"
# 
# INDEX="/bulk/jpoland/genome/wheat/Landmark/index/170831_Landmark_pseudomolecules_v1"
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/CDC-Landmark-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/CDC-Landmark-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/CDC-Landmark-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/CDC-Landmark-CENH3_R2.trim.fq"


###################################################################################################################################

## Mace ##

# NAME="mace"
# 
# INDEX="/bulk/jpoland/genome/wheat/Mace/index/181120_mace_pseudomolecule_v1.0"
# 
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/Mace-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/Mace-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/Mace-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/Mace-CENH3_R2.trim.fq"
# 

###################################################################################################################################

## Norin 61 ##

NAME="norin61"

##INDEX="/bulk/jpoland/genome/wheat/Norin61/index/180827_Norin61_pseudomolecule_v15"
INDEX="/bulk/jpoland/genome/wheat/Norin61/index/190307_Norin61_pseudomolecules_v1.1"

CHIP1="/bulk/jpoland/centromeres/CHIP002/Norin61-CENH3_R1.fq"
CHIP2="/bulk/jpoland/centromeres/CHIP002/Norin61-CENH3_R2.fq"

CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/Norin61-CENH3_R1.trim.fq"
CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/Norin61-CENH3_R2.trim.fq"


###################################################################################################################################

## CDC Stanley ##

# NAME="stanley"
# 
# INDEX="/bulk/jpoland/genome/wheat/Stanley/index/180902_Stanley_pseudomolecules_v1.2"
# 
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/CDC-Stanley-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/CDC-Stanley-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/CDC-Stanley-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/CDC-Stanley-CENH3_R2.trim.fq"
# 


###################################################################################################################################

## SY Mattis ##
# 
# NAME="mattis"
# 
# INDEX="/bulk/jpoland/genome/wheat/SY_Mattis/index/181016_SY_Mattis_pseudomolecule_v1"
# 
# CHIP1="/bulk/jpoland/centromeres/CHIP002/SY-Mattis-CENH3_R1.fq"
# CHIP2="/bulk/jpoland/centromeres/CHIP002/SY-Mattis-CENH3_R2.fq"
# 
# CHIP1TRIM="/bulk/jpoland/centromeres/CHIP002/SY-Mattis-CENH3_R1.trim.fq"
# CHIP2TRIM="/bulk/jpoland/centromeres/CHIP002/SY-Mattis-CENH3_R2.trim.fq"
# 



###################################################################################################################################




## 1R analysis ##

## NAME="larry"

## build index
# FASTA="/bulk/jpoland/genome/wheat-rye/index/wheat_161010_Chinese_Spring_v1.0-rye_Lo7_2018v1p1p1_1R_pseudomolecules.fasta"
# INDEX="/bulk/jpoland/genome/wheat-rye/index/wheat_161010_Chinese_Spring_v1.0-rye_Lo7_2018v1p1p1_1R_pseudomolecules"
 
# /homes/jpoland/tools/hisat2-2.1.0/hisat2-build -p 8 ${FASTA}  ${INDEX}

# R1="/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03708_HNC55CCXY_L4_1.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03708_HNC55CCXY_L8_1.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HNC55CCXY_L8_1.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HNYCGCCXY_L3_1.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HTGTWCCXY_L3_1.fq.gz"
# R2="/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03708_HNC55CCXY_L4_2.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03708_HNC55CCXY_L8_2.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HNC55CCXY_L8_2.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HNYCGCCXY_L3_2.fq.gz,/bulk/jpoland/genome/wheat/Larry/raw_data/PE/Larry_DDP03709_HTGTWCCXY_L3_2.fq.gz"

# SAM="/bulk/jpoland/centromeres/larry.sam"
# BAM="/bulk/jpoland/centromeres/larry.bam"
# BAMSORT="/bulk/jpoland/centromeres/larry.sort.bam"
# BAMFILTER="/bulk/jpoland/centromeres/larry.sort.filter.bam"
# 




###################################################################################################################################

# NAME="tam112"
# 
# INDEX="/bulk/jpoland/genome/wheat-rye/index/wheat_161010_Chinese_Spring_v1.0-rye_Lo7_2018v1p1p1_1R_pseudomolecules"
# 
# 
# R1="/bulk/jpoland/genome/wheat/TAM112/raw_data/PE/TAM112_DDP03710_HNC55CCXY_L6_1.fq.gz"
# R2="/bulk/jpoland/genome/wheat/TAM112/raw_data/PE/TAM112_DDP03710_HNC55CCXY_L6_2.fq.gz"
# 
# BULK="/bulk/jpoland/centromeres/"
# SAM="${BULK}${NAME}.sam"

## RUN HISAT2
# /homes/jpoland/tools/hisat2-2.1.0/hisat2 -p 8 -x ${INDEX} -1 ${R1} -2 ${R2} -S ${SAM} 






###################################################################################################################################




## trim adapters
# ~/.local/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 
# grep 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' /bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_SRR1686799/SRR1686799_2.fastq
# zcat /bulk/jpoland/genome/wheat/ChineseSpring/raw_data/ChipSeq/ChipSeqCENH3_KSU-Koo/KSU_AA-ChipSeqCENH3_NoIndex_L001_R1_001.fastq.gz | head -n1000000 | grep 'AGATCGGAAGAGCA'


# module load Java
# module load R
# module load Perl
# 
# sample_name=""
# R1=""
# R2=""
# threads=""
# 
# ###################################################################################################################################
# 
# java -jar /homes/sshrest1/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 8 ${CHIP1} ${CHIP2} \
# ${sample_name}_1P.fq ${sample_name}_1U.fq ${sample_name}_2P.fq ${sample_name}_2U.fq \
# ILLUMINACLIP:/homes/sshrest1/bin/Trimmomatic-0.38/adapters/illuminaClipping.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 &> TRIM_${sample_name}.log
# 
# java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36


## set names ##

BULK="/bulk/jpoland/centromeres/"
SAM="${BULK}${NAME}.sam"
BAM="${BULK}${NAME}.bam"
BAMSORT="${BULK}${NAME}.sort.bam"
BAMFILTER="${BULK}${NAME}.sort.filter.bam"

BINS="/homes/jpoland/centromeres/chr-100Kb-bin.bed"
## BINS="/homes/jpoland/centromeres/chr.wheat.1R-100Kb-bin.bed"


INTERSECT="/homes/jpoland/centromeres/${NAME}.intersect.bed"





## trim adapters ## --reads_to_process 100000000
# /homes/jpoland/tools/fastp/fastp --in1 ${CHIP1} --in2 ${CHIP2}  --out1 ${CHIP1TRIM} --out2 ${CHIP2TRIM} --thread 8 
# 
# mv fastp.json ${NAME}.fastp.json
# mv fastp.html ${NAME}.fastp.html

## RUN HISAT2
# /homes/jpoland/tools/hisat2-2.1.0/hisat2 -p 8 -x ${INDEX} -1 ${CHIP1TRIM} -2 ${CHIP2TRIM} -S ${SAM} 

## convert to BAM
# /homes/jpoland/tools/samtools/bin/samtools view -b -S -@ 8 -o ${BAM} ${SAM}

## sort
# /homes/jpoland/tools/samtools/bin/samtools sort -m 4G -@ 8 -o ${BAMSORT} ${BAM}

## filter
# /homes/jpoland/tools/samtools/bin/samtools view -q 30 -o ${BAMFILTER} ${BAMSORT}


## intersect bed
module load GLib
/homes/jpoland/tools/bedtools-2.29/bin/bedtools intersect -a ${BINS} -b ${BAMFILTER} -c -sorted > ${INTERSECT}



## rename the output files

echo 'nice job!'

# mv *${NAME}.o ${NAME}.out
# mv *${NAME}.e ${NAME}.error
# 





