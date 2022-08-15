# QVG

**Quick Viral genome Genotyper**

If you find this pipeline useful, please, do not forget to credit our work by citing this paper:

Varadi, A., Kaszab, E., Kardos, G., Prepost, E., Szarka, K., & Laczko, L. (2022). Rapid genotyping of viral samples using Illumina short-read sequencing data. *bioRxiv preprint*. https://doi.org/10.1101/2022.03.21.485184

Feedback is very welcome.

## Overview

This pipeline is designed to genotype targeted viral strains quickly. The input for the pipeline are the single or paired-end reads produced by Illumina sequencing experiments in `fastq.gz` format and the list of sample names to be analyzed, and a reference genome. Currently, the pipeline does not include any specific step to remove host contamination and therefore assumes that the target viral DNA is present at the highest frequency. After aligning the filtered reads to the preferably closely related reference genome, marking PCR artifacts and optical duplicates, and clipping high-depth alignments, two variant calls are performed. One of them assumes that the ploidy of samples is one and exports the sequence of the dominant genome found in each sample file, whereas the second assumes pooled sequencing and, after calculating the allele balance, gives insight into the heterogeneity (within-host variation) of samples. Optionally annotating the genomes or smoothing of sequencing depth is also possible; then, the genome annotation `.gff` file should also be provided.

## Installation and dependencies

The recommended installation of QVG is to clone this repository with the following command from GitHub: 

```
git clone https://github.com/laczkol/QVG.git 
```

After cloning, using the `conda` package manager, dependencies specified in ` qvg-env.yaml ` file can be easily installed by navigating to the copied directory and creating a `conda` environment for QVG:

```
cd ./QVG/
conda create --name qvg-env --file qvg-env.yaml 
```

Creating a conda environment ensures that dependencies can be found correctly and eliminates conflicts between different software versions; thus, this is the preferred way to install QVG. If Miniconda is not yet installed, please, visit [this site](https://docs.conda.io/en/latest/miniconda.html) to obtain it.
After that, the newly created environment with all dependencies installed can be activated by typing the following in the terminal:

````
conda activate qvg-env
````

This way, all dependencies except R and R Script will be installed and added to `$PATH`. R should be installed manually (see Dependencies).

To run QVG from the command line using a GNU/Linux operating system, the path to `QVG.sh` should be added to the user's $PATH variable. To do this, open the user's shell-specific configuration file (usually `~/.bashrc`) with a text editor (`nano` in this example).

```
nano ~/.bashrc
```

Add the following line to the end of it. Change the `path/to/script` part to the actual path to QVG.sh.

```
export PATH="path/to/QVG.sh:$PATH"
```

Save the file, then load the new `$PATH` to the current shell session:

```
source ~/.bashrc
```

For a temporary effect, you can `export PATH="path/to/QVG.sh:$PATH"` in your terminal, in which case the path of `QVG.sh` should be exported every time a new session is opened.

Dependencies are the following (the way to install them one by one is given in code blocks):

- [fastp](https://github.com/OpenGene/fastp) for filtering the raw reads

- [bwa](http://bio-bwa.sourceforge.net/) to align the short reads to a reference genome

- [samtools](http://www.htslib.org/) and [HTSlib](http://www.htslib.org/) for `sam` to `bam` conversion and to output alignment statistics.

- [sambamba](https://lomereiter.github.io/sambamba/) to mark duplicate alignments

- [freebayes](https://github.com/freebayes/freebayes) for variant calling

- [bcftools](https://samtools.github.io/bcftools/bcftools.html) to export variants and allele balance in tabular format

- vcf2fasta, vcfstats and vcffilter from the [vcflib](https://github.com/vcflib/vcflib) for file conversion, exporting vcf statistics and variant filtering

- [vcftools](http://vcftools.sourceforge.net/) to filter for missingness

- [bedtools](https://bedtools.readthedocs.io/en/latest/) to mask the genomic regions with no reads

- [liftoff](https://github.com/agshumate/Liftoff) and [minimap2](https://github.com/lh3/minimap2
) to annotate genomic regions

- [R](https://www.r-project.org/) and [Rscript](https://rdrr.io/r/utils/Rscript.html) for visualization. Installation of R is detailed its own website. Most GNU/Linux distributions can install `R` using its official repository.

- [GNU parallel](https://www.gnu.org/software/parallel/) to run tasks in parallel

- [bioawk](https://github.com/lh3/bioawk) for the basic processing of fasta sequences.

- [GNU Core Utilities](https://www.gnu.org/software/coreutils/) is practically the spine of the pipeline. It is used in most `bash` scripts and should be pre-installed on most GNU/Linux-based operating systems.

Software installed with conda is added automatically to the `$PATH` variable if `miniconda` is configured correctly.

All dependencies can be installed in the current environment by copying the following lines to the terminal:

```
conda install -y -c bioconda fastp
conda install -y -c bioconda bwa
conda install -y -c bioconda samtools=1.15.1
conda install -y -c bioconda sambamba=0.8.2
conda install -y -c bioconda freebayes
conda install -y -c bioconda bcftools
conda install -y -c bioconda vcflib
conda install -y -c bioconda vcftools
conda install -y -c bioconda bedtools
conda install -y -c bioconda liftoff minimap=2.17
conda install -y -c conda-forge parallel
conda install -y -c bioconda bioawk
```

NOTE: Copying the above code-block will install the dependencies one by one. Installing them in a single command with channels correctly specified the same effect should be achieved (`conda install -y -c bioconda -c conda-forge fastp bwa samtools=1.15.1 sambamba=0.8.2 freebayes bcftools vcflib vcftools bedtools liftoff minimap=2.17 parallel bioawk`). **Using the option `-y' will install the required tools in the current environment without the need to confirm the action. Please, make sure that installing the above dependencies does not interfere with other tools used on the same computer in the same environment or create a virtual environment for QVG as shown above.**

`R` and `Rscript` are not included in the provided `.yaml` file and should be installed manually. This is because even the newest R version installed with conda might have dependency issues at some systems (conflict of dependencies). Please ensure that `R` is installed correctly and added to your `$PATH` variable. A similar issue sometimes can be observed with samtools. Please ensure that typing `samtools` to your terminal does not throw any errors. In the provided `.yaml` file, `samtools 1.15.1` is included and is recommended.

After cloning the repository and installing the dependencies, it is advised to check if all dependencies can be found correctly and if the pipeline works as expected. The repository contains a small subset of SARS-CoV2 sequencing reads and a reference genome with the corresponding annotation. To check if the pipeline works using these test data, please, run `./run_test.sh` from the directory where the repository was cloned. If the output correctly tells when the run ended, QVG should work correctly on real data. If a line starting with `Run ended` is not output to the screen, the pipeline stopped somewhere during the analysis. Additionally, this test looks for some main output files of the pipeline, namely, the consensus genome sequence, the sites variable within-host, and the transferred genome annotation. If the test tells all these files were found, the installation of QVG is correct; otherwise, the availability of dependencies should be double-checked.

The help menu of the pipeline can be invoked by typing `QVG.sh -h` or `QVG.sh --help`, whereas the pipeline's version can be checked by typing `QVG.sh -v` or `QVG.sh --version`.

## Details

The pipeline can be parametrized from the command line. The two mandatory options to run the pipeline are the following:

```
 -r or --reference-genome
 	The reference genome sequence in .fasta format. The reference should contain only one contig.
 -samples-list or --samples-list
	A text file listing sample file basenames to be included in the analysis.
```

Please ensure that the file specified with `-samples-list` contains only the basenames of sample files and does not contain any headers, extra new-line characters, and empty lines. For example, this text file to include paired-end reads of five samples named S1, S2, S3, S4, and S5 would look like the following:

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

With all other parameters left by default, the pipeline can be invoked by typing the following line in a terminal (assuming that `QVG.sh` is added to the `$PATH`):

````bash
QVG.sh -r reference_genome.fasta -samples-list list_of_samples.txt
````

By default, the input and output directories are assumed to be the current directories. The following options can alter these parameters:

````
-s or --samples-directory
	The input directory containing the sample files listed in list_of_samples.txt.
-o or --output-directory
	The output directory to store all the output files of the pipeline.
````

QVG uses the GNU Parallel tool to speed up the analysis. The number of CPU threads used during the pipeline can be specified by `-np or number-of-processors` and the default value is one. For example, setting `-np 10` will use 10 CPU threads.

Other parameters that can be set and make the fine-tuning of the pipeline possible:

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
-trim_front1 or --trim-front1
	The number of bases to be trimmed from the beginning of R1 reads. [default = 10]
-trim_front2 or --trim-front2
	The number of bases to be trimmed from the beginning of R2 reads. [default = 10]
-trim_tail1 or --trim-tail1
	The number of bases to be trimmed from the end of R1 reads. [default = 10]
-trim_tail2 or --trim-tail2
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

The pipeline, by default, clips alignments at high-coverage regions that can be fine-tuned with the following options:

````
-cw or --clip-window
	The sliding window size used to assess read depth. [default = 100]
-cs or --clip-step
	The step size of sliding windows. [default = 10]
-hc or --high-coverage
	The mean read depth is assessed for each sample file. This value is used to multiply the mean read depth and define the read depth threshold of regions to be clipped. [default = 10]
````

Optionally, the read depth of alignments can be smoothed out along the reference genome (can be useful to decrease read depth bias). The following options can parametrize this feature: 

````
-sc or --smooth-coverage
	Defines if coverage smoothing should be done. [default = no]
-smoothw or --smooth-window
	The consecutive window size of random resampling. [default = 100]
-scount or --smooth-count
	Read count within the window for resampling. [default = 500]
````

Optionally, a `.gff` file can be specified to transfer the annotations of the reference. If a `.gff` file is supplied annotation transfer is automatically turned on.

````
-g or --gff-file
	The .gff file containing the gene annotations of the reference.
````

The screening of within-host variable sites can be turned off by the following option:

````
-pool or --pooled-sequencing
  Valid options are yes or no. Turns on or off the screening of within-host variability. If no such sites are expected it is advised to turn off this feature. [default yes]
````

A typical command to run `QVG.sh` including the annotation step, would look like the following:

````
QVG.sh -r reference_genome.fasta -samples-list list_of_samples.txt -s ./fastq_files -o ./output_files -annot yes -g reference_genome.gff3 -np <number_of_threads>
````

The default values generally work well with AmpliSeq data obtained by paired-end sequencing using Illumina Miseq. 

## Main output

The output directory should contain the following files:

`coverages.pdf`: This is the summary of statistics exported by `samtools`, including the number of positions and percent of the genome covered by sequencing reads,  the total number of reads, and the mean read depth. This plot is created from the table `coverages.tsv`.

`non_haploid_sites.txt`: The distribution of sites that have more than one probable allele assuming pooled sequencing. The first column represents the number of samples with a given ‘non-haploid’ site; the second column shows the position of such sites on the reference genome.

The dominant genome of samples is exported in multifasta format; these filenames contain the date and time of the run.

`samples_multifasta_YEAR_MONTH_DAY_TIME.fa`: The dominant genome of samples in multifasta format before masking regions with no reads.

`samples_multifasta_masked_YEAR_MONTH_DAY_TIME.fa`: The dominant genome of samples in multifasta format after masking regions with no reads.

Main output directories contain subdirectories with the basename of each sample file. These directories include all files created for the given sample during the run, including the quality-filtered reads, the short-read alignment with marked duplicates in `.bam` format, and statistics of sample files.

`BASENAME_filter.html` and `BASENAME_filter.json`: These contain the summary of the quality filtering as output by `fastp`. Filtered reads are stored in `BASENAME_R1.fq.gz` and `BASENAME_R2.fq.gz`.

`BASENAME_depth.png` and `BASENAME_depth.log.png`: These plots show the absolute value and log-transformed value of read depth along the reference genome. This is created from `BASENAME_depth.tsv`.

`BASENAME_coverage.tsv`, `BASENAME_coverage.txt`, `BASENAME_depth.tsv`, `BASENAME_flagstat.txt`, `BASENAME_idxstats.tsv` contain the output of the corresponding statistic as output by `samtools`.

`BASENAME_p1.vcf` and `BASENAME_p1_GT.tsv`: Genotype information are stored in these files in `vcf` and `tsv` format assuming a given ploidy (one by default).

`BASENAME_SNPdensity.snpden` and `BASENAME_SNPdensity.snpden.pdf`: These files show the SNP-density along the reference genome in tabular and in `pdf` format.

`BASENAME_vcfstats_p1.txt`: Assuming a given ploidy vcf statistics are stored in this file, which contains the total number of sites and the type of variants.

`BASENAME_pooled.vcf`, `BASENAME_pooled_GT.tsv` and `BASENAME_pooled_GT.pdf`: The variants observed when assuming pooled sequencing are stored in `vcf` and `tsv` format, whereas the frequency of allele-balance values are plotted and output in `pdf`.

`BASENAME_vcfstats_pooled.txt`: Assuming pooled sequencing vcf statistics are stored in this file, which contains the total number of sites and the type of variants.

The fasta file `BASENAME_REFERENCE-SEQUENCE-ID.fa` contains the sequence of the dominant genome found in the given sample before the masking of non-sequenced regions, and the masked genome sequence is stored in `BASENAME_masked.fasta`.
