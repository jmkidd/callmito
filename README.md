# callmito

Pipeline for calling canine mitochondiral variation from aligned Illumina data


# Summary

This pipeline calls mitochondrial variation from aligned Illumina data. The input is 
a BAM/CRAM of paired-end Illumina reads aligned to the UU_Cfam_GSD_1.0 assembly (as in 
[dogmap](https://github.com/jmkidd/dogmap)). The output is a reconstructed fasta of the 
sample mitochondrial genome and a VCF file of potential variants. 

Variant calling uses Mutect2 as described in [Laricchia et al.](https://genome.cshlp.org/content/early/2022/01/24/gr.276013.121.abstract)
and utilizes filters and assignment suggested by [Fregel et al.](https://pubmed.ncbi.nlm.nih.gov/25869968/).
Parameters were optimized based on comparison with mitocondrial genomes reported in recently
published genomes with PacBio and Illumina sequence data. 

# Conceptual overview of pipeline

## Step 1: Extract mitochodnrial reads
Read pairs with at least run read that align to chrM or to Nuclear Mitochondiral Segments 
(NUMTs) that are >300bp with > 95% identity are extracted.  NuMTs were identified by
Fabian Ramos-Almodovar 

## Step 2: Align to standard and rotated mitochodnrial reference genome
NC_002008.4 is used as the reference mitochondiral genome. To account for the cirular nature
of the mitochondria, extracted reads are aligned NC_002008.4 as well as to a version of 
NC_002008.4 that has been rotated by 8 kbp. Resulting BAM files are sorted and duplicate marked.

## Step 3: Coverage assessment
Coverage along the mitochondrial genome is assessed using CollectHsMetrics from GATK. Depth is
merged from the reference and rotated alignments, with the first and last 4kbp taken from
the rotated alignment.

If the mean coverage is greater than 5,000, the reference and rotated BAMs are downsampled
to a depth of 5,000 using DownsampleSam from GATK.

## Step 4: Call mitochondrial variation using Mutect2
Variants are called from the standard and rotated BAM using mitochondrial mode of Mutect2.
The following options are used: --mitochondria-mode --max-reads-per-alignment-start 75 --max-mnp-distance 0 --annotation StrandBiasBySample

Filters are applied to each of the resulting VCF files using GATK FilterMutectCalls --mitochondria-mode

## Step 5: Merge and filter VCF files to make germline calls
The reference and rotated VCF files are merged, with the  first and last 4kbp taken from
the rotated alignment. The resulting variants are filtered to remove sites where the most
frequent alternative allele fails the strand_bias filter or where the most frequenct alternative
allele has a allele fraction less than 0.5. 

## Step 6: Construct germline mitochondrial fasta 
A representation of the mitochondria sequence is created using bcftools consensus. Regions with
a coverage less than 100, or that overlap with NC_002008.4 positions 15990-15990 or 15512-15535
are masked to 'N'

## Step 7: Attempt haplogroup assignment
An attempt at haplogroup assignment based on the diagnostic changes reported in Fregel et al.
is reported. Note that these assignments should not be considered definitive due to ambiguities
and limited resolution. Assessment relative to clade defining haplogroup sequences is recommended.

# Software required

The following software and versions are used
```
bwa-mem
gatk version 4.2.5.0
samtools version >= 1.9
liftOver
bgzip
tabix
bcftools version >= 1.9.
```

# Usage
```
python callmito/process-sample.py \
--ref GENEOME-FASTA-USED-IN-ALIGNMENT \
--name SAMPLE-NAME \
--finaldir OUTPUT-DIRECTORY \
--cram SAMPLE-BAM-or-CRAM-FILE \
--coords COORDINTATES-TO-EXTRACT (callmito/numt-coords-to-extract.95_300.bed) \
--mitoFa MITO-REFERENCE (callmito/refs/NC_002008.4.fa) \
--mitoFaRotated ROTATED-MITO-REFERENCE (callmito/refs/NC_002008.4.rotate8k.fa) \
--chainfile LIFTOVER-CHAIN-FILE-FOR-ROTATED-MITO (callmito/refs/rotatedToOriginal.liftOver) \
--diagnosticTable  DIAGNOSITIC-SNPS-TABLE (callmito/fregel-haplogroups.txt)
```




