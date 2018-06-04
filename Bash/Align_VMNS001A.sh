#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Charles/CM_WGBS_WD_Mouse_Liver/
#SBATCH --mem=100000
#SBATCH --time=1-0:00
#SBATCH --partition=gc,gc64,gc128,gc256,gc512
#SBATCH --mail-type=END                     
#SBATCH --mail-user=cemordaunt@ucdavis.edu  

# Load Modules
PATH=$PATH:/share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/
module load bowtie/1.1.1
module load samtools/0.1.19
module load sratoolkit/2.4.2-3
module load bedtools2/2.25.0
export PYTHONPATH=/share/lasallelab/pysam/lib/python2.7/site-packages/

# Unzip and Filter VMNS001A
gunzip -c raw_files/VMNS_1A*fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > raw_files/VMNS001A_filtered.fq
gzip raw_files/VMNS001A_filtered.fq
perl /share/lasallelab/programs/perl_script/adapter_split.pl raw_files/VMNS001A_filtered.fq.gz raw_files/VMNS001A_noadap.fq.gz raw_files/VMNS001A_withadap.fq.gz
perl /share/lasallelab/programs/perl_script/adapter_trimmer.pl raw_files/VMNS001A_withadap.fq.gz raw_files/VMNS001A_trimmed.fq.gz 45 10

# Align VMNS001A to mm10
mkdir VMNS001A
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/mm10/mm10.fa -d /share/lasallelab/genomes/mm10/BSseek2_refgen/ -i raw_files/VMNS001A_noadap.fq.gz -o VMNS001A/VMNS001A_noadap.bam
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 80 -m 2 -f bam -g /share/lasallelab/genomes/mm10/mm10.fa -d /share/lasallelab/genomes/mm10/BSseek2_refgen/ -i raw_files/VMNS001A_trimmed.fq.gz -o VMNS001A/VMNS001A_trimmed.bam
samtools sort VMNS001A/VMNS001A_noadap.bam VMNS001A/VMNS001A_noadap_sorted
samtools sort VMNS001A/VMNS001A_trimmed.bam VMNS001A/VMNS001A_trimmed_sorted
samtools merge VMNS001A/VMNS001A.bam VMNS001A/VMNS001A_noadap_sorted.bam VMNS001A/VMNS001A_trimmed_sorted.bam
samtools view VMNS001A/VMNS001A.bam > VMNS001A/VMNS001A.sam

# Make PerMeth Files
mkdir VMNS001A/tmp
perl /share/lasallelab/programs/perl_script/SAMsorted_to_permeth.pl VMNS001A/VMNS001A.sam VMNS001A/tmp/PerMeth_VMNS001A VMNS001A mm10 CG combined 1
mkdir VMNS001A/PerMeth_VMNS001A
perl /share/lasallelab/programs/perl_script/gbcompliance.pl mm10 VMNS001A/tmp/PerMeth_VMNS001A_ VMNS001A/PerMeth_VMNS001A/PerMeth_VMNS001A_ VMNS001A VMNS001A
rm -r VMNS001A/tmp

# Make DSS files
perl /share/lasallelab/programs/perl_script/Permeth_to_DSSformat.pl mm10 VMNS001A/PerMeth_VMNS001A/PerMeth_VMNS001A_ DSS_files/VMNS001A_
