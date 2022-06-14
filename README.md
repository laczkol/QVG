# QVG

**Quick Viral genome Genotype**

If you find this pipeline useful, please, do not forget to credit our work by citing this paper:

Varadi, A., Kaszab, E., Kardos, G., Prepost, E., Szarka, K., & Laczko, L. (2022). Rapid genotyping of viral samples using Illumina short-read sequencing data. *bioRxiv preprint*. https://doi.org/10.1101/2022.03.21.485184

Feedback is very welcome.

## Overview

This pipeline is desgined to quickly genotype viral strains. The input for the pipeline are the reads produced by Illumina sequencing experiments in `fastq.gz` format and the list of sample names to be analyzed. After aligning the filtered reads to the preferably closely related reference genome, two variant calls are performed. One of them assumes that the ploidy of samples is one and export the sequence of the dominant genome found in each sample file, whereas the second assumes pooled sequencing and after calculating the allele balance gives insight into the heterogeneity of samples.

## Installation and dependencies

With dependencies installed and added to the `$PATH` variable, QVG can be installed by cloning this repository (`git clone https://github.com/laczkol/QVG.git `) and run from the command line using a GNU/Linux operating system.

Dependencies are the following:

- [fastp](https://github.com/OpenGene/fastp) for filtering the raw reads

- [bwa](http://bio-bwa.sourceforge.net/) to align the short reads to a reference genome

- [samtools](http://www.htslib.org/) and [HTSlib](http://www.htslib.org/) for `sam` to `bam` conversion and to output alignment statistics. Please, make sure that the samtools version used is > 1.10.

- [sambamba](https://lomereiter.github.io/sambamba/) to mark duplicate alignments

- [freebayes](https://github.com/freebayes/freebayes) for variant calling

- [bcftools](https://samtools.github.io/bcftools/bcftools.html) to export variants and allele balance in tabular format

- vcf2fasta, vcfstats and vcffilter from the [vcflib](https://github.com/vcflib/vcflib) for file conversion, exporting vcf statistics and variant filtering

- [vcftools](http://vcftools.sourceforge.net/) to filter for missingness

- [bedtools](https://bedtools.readthedocs.io/en/latest/) to mask the genomic regions with no reads

- [liftoff](https://github.com/agshumate/Liftoff) to annotate genomic regions

- [R](https://www.r-project.org/) and [Rscript](https://rdrr.io/r/utils/Rscript.html) for visualization

- [GNU parallel](https://www.gnu.org/software/parallel/) to run tasks in parallel

- [GNU Core Utilities](https://www.gnu.org/software/coreutils/) is practically the spine of the pipeline. It is used in the majority of `bash` scripts and should be preinstalled on most of the GNU/Linux-based operating systems.

  After cloning this repository, using the `conda` package manager dependencies can be installed by navigating to the copied directory and typing this line into the terminal:

  ```
  conda create --name qvg-env --file qvg-env.yaml 
  ```

After that, the newly created environment that has all dependencies installed can be activated by typing the following in the terminal:

````
conda activate qvg-env
````

R and Rscript are not included in the provided `.yaml` file and should be installed manually. The reason for this is that even the newest version of R installed with conda might have dependency issues at some systems. Please, make sure that R is installed correctly and is added to your $PATH variable. A similar issue sometimes can be observed with samtools. Please, make sure that typing `samtools` to your terminal does not throw any errors. In the provided `.yaml` file samtools 1.15.1 is included and is recommended to run the pipeline, but any version of this software above v1.10 should work correctly.

## Details

The pipeline can parametrized from the command line. The two mandatory options to run the pipeline are the following:

```
 -r or --reference-genome
 	The file that contains the reference genome that should be used for both aligning the reads haplotype calling. The reference should contain only one contig.
 -samples-list or --samples-list
	Text file that list of sample file basenames to be included in the analysis
```

Please, make sure that file specified with `-samples-list` contains only the basenames of sample files and does not contain any headers, extra new-line characters and empty lines. For example, this text file to include paired-end reads of five samples named S1, S2, S3, S4, S5 would look like the following:

````
S1
S2
S3
S4
S5
````

By the naming convention applied, the pipeline would look for the following `fastq.gz` files:

````
S1_*R1*.fastq.gz	S1_*R2*.fastq.gz
S2_*R1*.fastq.gz	S2_*R2*.fastq.gz
S3_*R1*.fastq.gz	S3_*R2*.fastq.gz
S4_*R1*.fastq.gz	S4_*R2*.fastq.gz
S5_*R1*.fastq.gz	S5_*R2*.fastq.gz
````

The asterisk (*) denotes a wild-card and can be any character.

With all other parameters left default, the pipeline can be invoked by typing the following line in a terminal (assuming that `QVG.sh` is added to the `$PATH`):

````bash
QVG.sh -r reference_genome.fasta -samples-list list_of_samples.txt
````

By default, both the input and output directories are assumed to be the current directory. These paramteres can be altered by the following paramters:

````
-s or --samples-directory
	The input directory that contains the sample files listed in list_of_samples.txt.
-o or --output-directory
	The output directory that will store all the output files of the pipeline.
````

QVG uses the GNU Parallel tool to speed up the analysis. The number of CPU threads to be used during the whole pipeline can be specified by `-np or number-of-processors` and the default value is. For example, setting `-np 10` will use 10 CPU threads.

Other parameters that can be set and makes the fine-tuning of the pipeline possible:

````
-bwa_k or --min-seed-length
	The minimum seed length parameter of bwa. Used during short-read alignment. [default = 19]
-bwa_A or --matching-score
	The matching score parameter of bwa. Used during short-read alignment. [default = 1]
-bwa_B or --mismatch-penalty
	The mismatch penalty score of bwa. Used during short-read alignment. [default = 4]
-bwa_O or --gap-open-penalty
	The gap opening penalty score of bwa. Used during short-read alignment. [default = 6]
-type or --sequencing-type
	The pipeline can use both singe-end and paired-end reads. The default is paired-end (PE). Setting this parameter to SE will look for single-end read files.
-trim_front1 or --trim_front1
	The number of bases to be trimmed from the beginning of R1 reads. [default = 10]
-trim_front2 or --trim_front2
	The number of bases to be trimmed from the beginning of R2 reads. [default = 10]
-trim_tail1 or --trim_tail1
	The number of bases to be trimmed from the end of R1 reads. [default = 10]
-trim_tail2 or --trim_tail2
	The number of bases to be trimmed from the end of R2 reads. [default = 10]
-minlen or --min-read-lenght
	The minimum length of reads to be included after quality filtering. [default = 30]
-p or --fb-ploidy
	Ploidy assumed for variant calling. [default = 1 for the first variant calling]
-Q or --fb-mismatch-base-quality-threshold
	The quality of the mismatched base to be included in the variant calling. [default = 30]
-m or --fb-min-mapping-quality
	The minimum mapping quality for a read to be included in the variant calling. [default = 30]
-mc or --fb-min-coverage
	Minimum coverage of a variant to be considered. [default = 5]
-n or --fb-best-alleles
	Number of most probable alleles to be considered for variant calling. [default = 5]
-F or --fb-min-alternate-fraction
	Minimal supporting fraction of reads for an alternate allele. [default = 0.2]
-mvq or --min-variant-quality
	Minimum quality of a variant to be kept. [default = 10]
-mqa or --min-qual-ao
	Minimum ratio of variant quality and observation count to be kept. [default = 10]
-sw or --snp-window
	Size of the sliding window to check for SNP-density. [default = 1000]
-mincov or --minimum-coverage
	The percent of the reference genome that should be covered to include a sample file in the analysis. [default = 95]
````
The pipeline by default clips alignments at high-coverage regions, that can be fine-tuned with the following options:
````
-cw or --clip-window
	The sliding window size used to assess read depth [default = 100]
-cs or --clip-step
	The step size of sliding windows [default = 10]
-hc or --high-coverage
	The mean read depth is assessed for each sample files. This value is used to multiply the mean read depth and define the read depth threshold of regions to be clipped [default = 10]
````
Optionally, the read depth of alignments can be smoothed out along the reference genome (can be useful to decrease read depth bias). This feature can be parametrized by the following options: 
````
-sc or --smooth-coverage
	Defines if coverage smoothing should be done [default = no]
-smoothw or --smooth-window
	The consecutive window size of random resampling [default = 100]
-scount or --smooth-count
	Read count within the window for resampling [default = 500]
````
Optionally, a gff file can be specified to transfer the annotations of the reference.

````
-annot or --annotate
	Valid options are yes or no. Specifies if annotation transfer should be carried out.
-g or --gff-file
	The gff file that contains the gene annotations of the reference.
````

These default values generally work well with AmpliSeq data obtained by paired-end sequencing using Illumina Miseq. 

## Main output

The output directory should containt the following files:

`coverages.pdf`: This is the summary of statistics exported by `samtools`, including the number of positions and percent of genome covered by sequencing reads,  the total number of reads and the mean read depth. This plot is created from the table `coverages.tsv`.

`non_haploid_sites.txt`: The distribution of sites that have more than one probable allels assuming pooled sequencing. The first column represents the number of samples a given ‘non-haploid’ site occurs and the second column shows the position of such sites on the reference genome.

The dominant genome of samples are exported in multifasta format and these filenames contain the date and time of the run.

`samples_multifasta_YEAR_MONTH_DAY_TIME.fa`: The dominant genome of samples in multifasta format before masking regions with no reads.

`samples_multifasta_masked_YEAR_MONTH_DAY_TIME.fa`: The dominant genome of samples in multifasta format after masking regions with no reads.

Main output directories containt subdirectories with the basename of each sample file. These directories include all files created for the given sample during the run, including the quality filtered reads , the short-read alignment with marked duplicates in `.bam` format, and statistics of sample files.

`BASENAME_filter.html` and `BASENAME_filter.json`: These contain the summary of the quality filtering as output by `fastp`. Filtered reads are stored in `BASENAME_R1.fq.gz` and `BASENAME_R2.fq.gz`.

`BASENAME_depth.png` and `BASENAME_depth.log.png`: These plots show the absolute value and log-transformed value of read depth along the reference genome. This is created from `BASENAME_depth.tsv`.

`BASENAME_coverage.tsv`, `BASENAME_coverage.txt`, `BASENAME_depth.tsv`, `BASENAME_flagstat.txt`, `BASENAME_idxstats.tsv` contain the output of the corresponding statistic as output by `samtools`.

`BASENAME_p1.vcf` and `BASENAME_p1_GT.tsv`: Genotype information are stored in these files in `vcf` and `tsv` format assuming a given ploidy (one by default).

`BASENAME_SNPdensity.snpden` and `BASENAME_SNPdensity.snpden.pdf`: These files show the SNP-density along the reference genome in tabular and in `pdf` format.

`BASENAME_vcfstats_p1.txt`: Assuming a given ploidy vcf statistics are stored in this file, which contain the total number of sites and the type of variants.

`BASENAME_pooled.vcf`, `BASENAME_pooled_GT.tsv` and `BASENAME_pooled_GT.pdf`: The variants observed when assuming pooled sequencing are stored in `vcf` and `tsv` format, whereas the frequency of allele-balance values are plotted and output in `pdf`.

`BASENAME_vcfstats_pooled.txt`: Assuming pooled sequencing vcf statistics are stored in this file, which contain the total number of sites and the type of variants.

The fasta file `BASENAME_REFERENCE-SEQUENCE-ID.fa` contains the sequence of the dominant genome found in the given sample before the masking of non-sequenced regions and the masked genome sequence is stored in `BASENAME_masked.fasta`.
