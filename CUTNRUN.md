CUT&RUN with Nextflow
================
Kathryn Lande
2026-05-07

### Please note: this was written on 7 May 2026 in reference to an incompatibility between the newest versions of nextflow (v25+) and the current version of nf-core/cutandrun which is on v3.2.2 as of today. The workaround in this document may no longer be necessary if nf-core/cutandrun has updated. Please check the nf-core/cutandrun [release history](https://nf-co.re/cutandrun/releases_stats/).

## This document is intended for use by IGC users.

# Run nf-core CUT&RUN

nf-core/cutandrun exclusively works on barrel at this time. To initiate
a run, you need to grab an older version of nextflow.

``` bash
# Grabbing a specific nextflow version:
curl -fsSL https://github.com/nextflow-io/nextflow/releases/download/v23.10.1/nextflow > nextflow_23.10.1
chmod +x nextflow_23.10.1

conda activate nf-core
```

Create a sample sheet, for cut and run it would look like this:

| group | replicate | fastq_1              | fastq_2              | control |
|-------|-----------|----------------------|----------------------|---------|
| Trt   | 1         | /path/S1_R1.fastq.gz | /path/S1_R2.fastq.gz | IgG     |
| Trt   | 2         | /path/S2_R1.fastq.gz | /path/S2_R2.fastq.gz | IgG     |
| Ctl   | 1         | /path/S3_R1.fastq.gz | /path/S3_R2.fastq.gz | IgG     |
| Ctl   | 2         | /path/S4_R1.fastq.gz | /path/S4_R2.fastq.gz | IgG     |
| IgG   | 1         | /path/S5_R1.fastq.gz | /path/S5_R2.fastq.gz |         |
| IgG   | 2         | /path/S6_R1.fastq.gz | /path/S6_R2.fastq.gz |         |

Run the pipeline on barrel with standard settings:

``` bash
./nextflow_23.10.1 run nf-core/cutandrun \
  --input SampleSheet.csv \
  --genome mm10 \
  -profile singularity \
  --blacklist /vast/igc/analyses/kat/Kat_Files/mm10_blacklist.bed \
  --gene_bed /vast/igc/analyses/kat/Kat_Files/mm10_genes.bed \
  --skip_fastqc true \
  --skip_igv true \
  --normalisation_mode CPM \
  --peakcaller macs2
```

## Get FRIPs

To call FRiPs manually, you can use
*/vast/igc/analyses/kat/Kat_Files/readFracs.sh*, which calls FRiPs
against sample-specific peaks, consensus peaks, and mm10 promoters for
all samples output by the nf-core/cutandrun pipeline. The full script
looks like this:

``` bash
#!/bin/bash

# initiate files
echo -e "Sample\tTotalReads\tSoloFRiP\tConsensusFRiP\tFracInPromoter" > ReadFractions.txt

# if using nextflow's Cut&Run pipeline, changing these paths should suffice. You may have to change the sample name regex otherwise.
bampath="results/02_alignment/bowtie2/target/markdup"
peakpath="results/03_peak_calling/04_called_peaks/macs2"
peakpathcon="results/03_peak_calling/05_consensus_peaks"

# individual peaks
for bam in $bampath/*.markdup.sorted.bam
do
    sample=${bam//".target.markdup.sorted.bam"}
    sample=${sample//$bampath}
    sample=${sample//"/"}
    
    #conSamp=${sample//"_R2"}
    #conSamp=${conSamp//"_R1"}
    conSamp=$(echo "$sample" | sed 's/_R[0-9]*//')
    echo "Group: " $conSamp " | Sample: " $sample
    
    echo "Getting input peaks..."
    scp $peakpathcon"/"$conSamp".macs2.consensus.peak_counts.bed" peaksCon.bed
    scp $peakpath"/"$sample".macs2_peaks.narrowPeak" peaks.bed
    
    # total number of peaks
    echo "Counting total reads..."
    tot=`samtools view -c $bam`
    
    # total number of reads in peaks
    echo "Converting bam to bed..."
    bedtools bamtobed -i $bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > bambed_tmp.bed
    
    echo "Counting solo overlapping reads..."
    inPeak=`bedtools sort -i peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a bambed_tmp.bed -b stdin | wc -l`
    frip=$(echo "scale=5 ; $inPeak / $tot" | bc)
    
    echo "Counting consensus overlapping reads..."
    inPeak=`bedtools sort -i peaksCon.bed | bedtools merge -i stdin | bedtools intersect -u -a bambed_tmp.bed -b stdin | wc -l`
    fripCon=$(echo "scale=5 ; $inPeak / $tot" | bc)
    
    echo "Counting mm10 TSS overlapping reads..."
    inPeak=`bedtools sort -i /vast/igc/analyses/kat/Kat_Files/mm10_promoters.bed | bedtools merge -i stdin | bedtools intersect -u -a bambed_tmp.bed -b stdin | wc -l`
    fritss=$(echo "scale=5 ; $inPeak / $tot" | bc)
    
    echo "Removing temporary files..."
    rm bambed_tmp.bed peaks.bed peaksCon.bed
    echo -e "\n"$sample"\t"$tot"\t"$frip"\t"$fripCon"\t"$fritss >> ReadFractions.txt
done
```

This script can be executed from the same directory where you executed
the nf-core/cutandrun pipeline (the directory containing the results
folder). **Note: You must be in the base environment for this to run.**

``` bash
bash /vast/igc/analyses/kat/Kat_Files/readFracs.sh
```
