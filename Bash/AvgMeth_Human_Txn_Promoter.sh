#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Charles/DK_WGBS_HumanLiver/
#SBATCH --mem=24000
#SBATCH -n 1
#SBATCH --time=0-4:00
#SBATCH --partition=gc,gc64,gc128,gc256,gc512
#SBATCH --mail-type=END                     
#SBATCH --mail-user=cemordaunt@ucdavis.edu 

perl /share/lasallelab/programs/perl_script/AvgMeth.2col.pl WD_Human_Thioredoxin_Promoter_AvgMeth2col.txt Human_Thioredoxin_Gene_Promoters.bed 3 1 1 VMDK1C/PerMeth_VMDK1C/PerMeth_VMDK1C_ VMDK001C VMDK1D/PerMeth_VMDK1D/PerMeth_VMDK1D_ VMDK001D VMDK4C/PerMeth_VMDK4C/PerMeth_VMDK4C_ VMDK004C VMDK5C/PerMeth_VMDK5C/PerMeth_VMDK5C_ VMDK005C VMDK6B/PerMeth_VMDK6B/PerMeth_VMDK6B_ VMDK006B VMDK6C/PerMeth_VMDK6C/PerMeth_VMDK6C_ VMDK006C VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_ VMDK001A VMDK1B/PerMeth_VMDK1B/PerMeth_VMDK1B_ VMDK001B VMDK2A/PerMeth_VMDK2A/PerMeth_VMDK2A_ VMDK002A VMDK2B/PerMeth_VMDK2B/PerMeth_VMDK2B_ VMDK002B VMDK3A/PerMeth_VMDK3A/PerMeth_VMDK3A_ VMDK003A VMDK3B/PerMeth_VMDK3B/PerMeth_VMDK3B_ VMDK003B VMDK4A/PerMeth_VMDK4A/PerMeth_VMDK4A_ VMDK004A VMDK4B/PerMeth_VMDK4B/PerMeth_VMDK4B_ VMDK004B VMDK5A/PerMeth_VMDK5A/PerMeth_VMDK5A_ VMDK005A VMDK5B/PerMeth_VMDK5B/PerMeth_VMDK5B_ VMDK005B VMDK2C/PerMeth_VMDK2C/PerMeth_VMDK2C_ VMDK002C VMDK2D/PerMeth_VMDK2D/PerMeth_VMDK2D_ VMDK002D VMDK3C/PerMeth_VMDK3C/PerMeth_VMDK3C_ VMDK003C VMDK3D/PerMeth_VMDK3D/PerMeth_VMDK3D_ VMDK003D VMDK6A/PerMeth_VMDK6A/PerMeth_VMDK6A_ VMDK006A 