## Introduction

**nf-core/rnaseq** pipeline is used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## The nextflow tower public hosted server is used to run this pipeline.

To run the pipeline, go to tower.nf.

## Files prepared to run the pipeline.

1. Transfer all the fastq.gz files to aws s3 use aws s3 cp command line.
2. generate sample sheet as describe in the [nf-core rnaseq usage document](https://nf-co.re/rnaseq/3.3/usage). In the samplesheet, you provde the sample name, the path to fastq file(s), and the strandness information. transfer this file to aws s3 bucket.
3. Now you are all set to run the pipeline at nf-tower: tower.nf
4. In the nf-tower, you will need to inut the following. Most other parameters have defaults, but you can change.

- inut: the path to sample sheet
- outdir: the results directory.  **don't forget the last slash /**
- genome.

```
{
    "input":"s3://dicerna-etl/test/monkey/sample_sheet_full.csv",
    "outdir":"s3://dicerna-etl/test/monkey_dicerna/results_star/",  
    "genome":"mf5_ens_v102",
    "aligner":"hisat2"
}
```
Those customized genomes were hosted in aws s3 for dicerna use. We can upgrade as we want. The GRCm38_v102, mf6_ens_v104, ... can be used in the genome field.

 ```
      'human_ensembl_v104_2021_March' {
        fasta = 's3://dicerna-sysbio/genomes/human_ensembl_v014_2021_March/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
        gtf  = 's3://dicerna-sysbio/genomes/human_ensembl_v014_2021_March/Homo_sapiens.GRCh38.104.gtf.gz'
      }
     'human_refseq_release_109_2021_May' {
        fasta = 's3://dicerna-sysbio/genomes/human_refseq_release_109_2021_May/GCF_000001405.39_GRCh38.p13_genomic.fna.gz'
        gtf  = 's3://dicerna-sysbio/genomes/human_refseq_release_109_2021_May/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz'
      }
      // mouse GRCm38, ensembl v102
      'GRCm38_v102' {
        fasta = 's3://dicerna-sysbio/genomes/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
        gtf  = 's3://dicerna-sysbio/genomes/GRCm38/Mus_musculus.GRCm38.102.gtf'
      }
      // monkey, mf6, ensembl v104
      'mf6_ens_v104' {
        fasta = 's3://dicerna-etl/genomes/Macaca_fascicularis_6.0/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz'
        gtf = 's3://dicerna-etl/genomes/Macaca_fascicularis_6.0/Macaca_fascicularis.Macaca_fascicularis_6.0.104.gtf.gz'
      }
      // monkey, mf5, ensemb, v102
      'mf5_ens_v102' {
        fasta = 's3://dicerna-etl/genomes/Macaca_fascicularis_5.0_ensembl_release_102/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa.gz'
        gtf = 's3://dicerna-etl/genomes/Macaca_fascicularis_5.0_ensembl_release_102/Macaca_fascicularis.Macaca_fascicularis_5.0.102.gtf.gz'
      }
      // monkey, mf5, refseq release 101
      'mf5_refseq_r101' {
        fasta = 's3://dicerna-etl/genomes/Macaca_fascicularis_5.0_refseq_release_101/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz'
        gtf  = 's3://dicerna-etl/genomes/Macaca_fascicularis_5.0_refseq_release_101/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fixed.gtf.gz'
      }
```
 
Again, this is the check list for the important things to check:
```
Input/output options
  input          : s3://dicerna-etl/test/monkey/sample_sheet_full.csv
  outdir         : s3://dicerna-etl/test/monkey_dicerna/results_star/
Reference genome options
  genome         : mf5_ens_v102
  aligner        : "star_salmon"  or "hisat2"
```
  
## Pipeline summary

The SRA download functionality has been removed from the pipeline (`>=3.2`) and ported to an independent workflow called [nf-core/fetchngs](https://nf-co.re/fetchngs). You can provide `--nf_core_pipeline rnaseq` when running nf-core/fetchngs to download and auto-create a samplesheet containing publicly available samples that can be accepted directly as input by this pipeline.

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
4. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
5. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
6. Choice of multiple alignment and quantification routes:
    1. [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
    2. [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
    3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
7. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
8. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
9. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
10. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
11. Create bigWig coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
12. Extensive quality control:
    1. [`RSeQC`](http://rseqc.sourceforge.net/)
    2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
    3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
    4. [`Preseq`](http://smithlabresearch.org/software/preseq/)
    5. [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
13. Pseudo-alignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/); *optional*)
14. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

> * **NB:** Quantification isn't performed if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments. However, you can use this route if you have a preference for the alignment, QC and other types of downstream analysis compatible with the output of HISAT2.
> * **NB:** The `--aligner star_rsem` option will require STAR indices built from version 2.7.6a or later. However, in order to support legacy usage of genomes hosted on AWS iGenomes the `--aligner star_salmon` option requires indices built with STAR 2.6.1d or earlier. Please refer to this [issue](https://github.com/nf-core/rnaseq/issues/498) for further details.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_. Note: This pipeline does not currently support running with Conda on macOS if the `--remove_ribo_rna` parameter is used because the latest version of the SortMeRNA package is not available for this platform.

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/rnaseq -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run nf-core/rnaseq \
        --input samplesheet.csv \
        --genome GRCh37 \
        -profile <docker/singularity/podman/conda/institute>
    ```

    * An executable Python script called [`fastq_dir_to_samplesheet.py`](https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py) has been provided if you would like to auto-create an input samplesheet based on a directory containing FastQ files **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```console
        wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py
        ./fastq_dir_to_samplesheet.py <FASTQ_DIR> samplesheet.csv --strandedness reverse
        ```

## Documentation

The nf-core/rnaseq pipeline comes with documentation about the pipeline [usage](https://nf-co.re/rnaseq/usage), [parameters](https://nf-co.re/rnaseq/parameters) and [output](https://nf-co.re/rnaseq/output).

## Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

The pipeline was re-written in Nextflow DSL2 by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

Many thanks to other who have helped out along the way too, including (but not limited to):
[@Galithil](https://github.com/Galithil),
[@pditommaso](https://github.com/pditommaso),
[@orzechoj](https://github.com/orzechoj),
[@apeltzer](https://github.com/apeltzer),
[@colindaven](https://github.com/colindaven),
[@lpantano](https://github.com/lpantano),
[@olgabot](https://github.com/olgabot),
[@jburos](https://github.com/jburos).
