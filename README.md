# Benchmarking experiments.
The code base is adapted from the benchmarking code of Pangenie published  in Nature Genetics doi: [10.1038/s41588-022-01043-w](https://www.nature.com/articles/s41588-022-01043-w)

# Prerequisites

## Tools to be installed:
1. pangenie (https://github.com/eblerjana/pangenie)
2. Graphtyper(https://github.com/DecodeGenetics/graphtyper)
3. paragraph (https://github.com/Illumina/paragraph)
4. Gatk (https://gatk.broadinstitute.org/hc/en-us)

## Input data
1. Reference Genome(GRCh38/hg38) https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
2. Download Haploid resolved Assemblies:
```
# Download HG002
wget  ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/NA24385-denovo-H1.fa.gz
wget  ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/NA24385-denovo-H2.fa.gz

#Download  HG00731
wget  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200417_Marschall-Eichler_NBT_hap-assm/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta
wget  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200417_Marschall-Eichler_NBT_hap-assm/HG00731_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta

```
2. Call Phased Variants from haploid resolved assemblies using this Snakemake workflow to call phased variants: https://bitbucket.org/jana_ebler/vcf-merging/src/master/

3. Download Short read sample for HG00731.
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/004/ERR3241754/ERR3241754_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/004/ERR3241754/ERR3241754_2.fastq.gz
```

4. Subsample the Short reads to 5,10,20, and 30x coverage, and map all subsamples to the reference genome using bwa.

5. Build CCDG using all the subsamples using snakemake script in build_test_CCDG/
    1. edit config.json to specify the outputfolder and kmersize(we used 31)
    2. Enter the pathes of the subsampled fastq files in subsample_table.csv
    3. run 'snakemake -j8 --use-conda'


# Run The Benchmarking workflow
The benchmarking snakemake workflow is in benchmark subfolder. We first need to edit config.json to configure the pathes of the input and output data. run the workflow using 'snakemake -j8 --use-conda'


