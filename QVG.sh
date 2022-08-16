#!/bin/bash
set -euo pipefail -euo nounset

hm="QVG (Quick Viral genome Genotyper) v.0.8.1\n\n
Required input files:\n
 -r or --reference-genome\n\t				The reference genome sequence in .fasta format. The reference should contain only one contig.\n\n
 -samples-list or --samples-list\n\t		A text file listing sample file basenames to be included in the analysis.\n\n

Optional parameters:\n
-s or --samples-directory\n\t				The input directory containing the sample files listed in list_of_samples.txt.\n\n
-o or --output-directory\n\t				The output directory to store all the output files of the pipeline.\n\n
-bwa_k or --min-seed-length\n\t				The minimum seed length parameter of bwa. Used during short-read alignment. [default = 19]\n\n
-bwa_A or --matching-score\n\t				The matching score parameter of bwa. Used during short-read alignment. [default = 1]\n\n
-bwa_B or --mismatch-penalty\n\t			The mismatch penalty score of bwa. Used during short-read alignment. [default = 4]\n\n
-bwa_O or --gap-open-penalty\n\t			The gap opening penalty score of bwa. Used during short-read alignment. [default = 6]\n\n
-type or --sequencing-type\n\t				The pipeline can use both singe-end and paired-end reads. The default is paired-end (PE). Setting this parameter to SE will look for single-end read files.\n\n
-trim_front1 or --trim-front1\n\t			The number of bases to be trimmed from the beginning of R1 reads. [default = 10]\n\n
-trim_front2 or --trim-front2\n\t			The number of bases to be trimmed from the beginning of R2 reads. [default = 10]\n\n
-trim_tail1 or --trim-tail1\n\t				The number of bases to be trimmed from the end of R1 reads. [default = 10]\n\n
-trim_tail2 or --trim-tail2\n\t				The number of bases to be trimmed from the end of R2 reads. [default = 10]\n\n
-minlen or --min-read-length\n\t			The minimum length of reads to be included after quality filtering. [default = 30]\n\n
-p or --fb-ploidy\n\t						Ploidy assumed for variant calling. [default = 1 for the first variant calling]\n\n
-Q or --fb-mismatch-base-quality-threshold\n\t	The quality of the mismatched base to be included in the variant calling. [default = 30]\n\n
-m or --fb-min-mapping-quality\n\t			The minimum mapping quality for a read to be included in the variant calling. [default = 30]\n\n
-mc or --fb-min-coverage\n\t				Minimum coverage of a variant to be considered. [default = 5]\n\n
-n or --fb-best-alleles\n\t					Number of most probable alleles to be considered for variant calling. [default = 5]\n\n
-F or --fb-min-alternate-fraction\n\t		Minimal supporting fraction of reads for an alternate allele. [default = 0.2]\n\n
-mvq or --min-variant-quality\n\t			Minimum quality of a variant to be kept. [default = 10]\n\n
-mqa or --min-qual-ao\n\t					Minimum ratio of variant quality and observation count to be kept. [default = 10]\n\n
-sw or --snp-window\n\t						Size of the sliding window to check for SNP-density. [default = 1000]\n\n
-mincov or --minimum-coverage\n\t			The percent of the reference genome that should be covered to include a sample file in the analysis. [default = 95]\n\n
-cw or --clip-window\n\t					The sliding window size used to assess read depth. [default = 100]\n\n
-cs or --clip-step\n\t						The step size of sliding windows. [default = 10]\n\n
-hc or --high-coverage\n\t					The mean read depth is assessed for each sample file. This value is used to multiply the mean read depth and define the read depth threshold of regions to be clipped. [default = 10]\n\n
-sc or --smooth-coverage\n\t				Defines if coverage smoothing should be done. [default = no]\n\n
-smoothw or --smooth-window\n\t				The consecutive window size of random resampling. [default = 100]\n\n
-scount or --smooth-count\n\t				Read count within the window for resampling. [default = 500]\n\n
-g or --gff-file\n\t						The .gff file containing the gene annotations of the reference. If provided the annotation transfer is automatically turned on.\n\n
-pool or --pooled-sequencing\n\t			Valid options are yes or no. Turns on or off the screening of within-host variability. If no such sites are expected it is advised to turn off this feature. [default yes]\n\n\n

-h or --help\n\t							This help menu.\n\n
-v or --version\n\t							Version information of QVG.\n\n\n

Example:\n
QVG.sh -r reference_genome.fasta -samples-list list_of_samples.txt -s ./fastq_files -o ./output_files -annot yes -g reference_genome.gff3 -np <number_of_threads>\n
"

cdir=$(pwd)
indir=$(pwd)
outdir=$(pwd)
ref_db=" "

np=1
annot="no"
gff=" "
slist=" "
call="no"
type="PE"
bwa_k=19
bwa_A=1
bwa_B=4
bwa_O=6
ploidy=1
minqual=30
mismatchqual=30
mapq=30
fb_min_cov=5
nbest=5
altfrac=0.2
tf1=10
tf2=10
tt1=1
tt2=1
ml=30
min_var_qual=10
min_qual_ao=10
snpwindow=1000
mincov=95
clipwindow=100
clipstep=10
hcm=10
smooth_coverage="no"
smooth_window=100
smooth_count=500
pool="yes"

while [[ "$#" -gt 0 ]];
	do
	case $1 in
		-s|--samples-directory)
			indir="$2"
			shift
			;;
		-o|--output-directory)
			outdir="$2" #output directory, USE ABSOLUTE PATH
			shift
			;;
		-r|--reference-genome) #reference genome in fasta, MANDATORY
			ref_db=$2
			shift
			;;
		-np|--number-of-processors) #number of processors
			np=$2
			shift
			;;
		-bwa_k|--min-seed-length) #minimum seed length of bwa
			bwa_k=$2
			shift
			;;
		-bwa_A|--matching-score) #matching score for bwa
			bwa_A=$2
			shift
			;;
		-bwa_B|--mismatch-penalty) #mismatch-penalty for bwa
			bwa_B=$2
			shift
			;;
		-bwa_O|--gap-open-penalty) #gap-open penalty for bwa
			bwa_O=$2
			shift
			;;
		-samples-list|--samples-list) #samples list MANDATORY
			slist=$2
			shift
			;;
		-type|--sequencing-type) #SE or PE,default is PE
			type=$2
			shift
			;;
		-trim_front1|--trim-front1)
			tf1=$2
			shift
			;;
		-trim_front2|--trim-front2)
			tf2=$2
			shift
			;;
		-trim_tail1|--trim-tail1)
			tt1=$2
			shift
			;;
		-trim_tail2|--trim-tail2)
			tt2=$2
			shift
			;;
		-minlen|--min-read-length)
			ml=$2
			shift
			;;
		-cw|--clip-window)
			clipwindow=$2
			shift
			;;
		-cs|--clip-step)
			clipstep=$2
			shift
			;;
		-hc|--high-coverage)
			hcm=$2
			shift
			;;
		-sc|--smooth-coverage)
			smooth_coverage=$2
			shift
			;;
		-smoothw|--smooth-window)
			smooth_window=$2
			shift
			;;
		-scount|--smooth-count)
			smooth_count=$2
			shift
			;;
		-p|--fb-ploidy)
			ploidy=$2
			shift
			;;
		-q|--fb-min-base-quality)
			minqual=$2
			shift
			;;
		-Q|--fb-mismatch-base-quality-threshold)
			mismatchqual=$2
			shift
			;;
		-m|--fb-min-mapping-quality)
			mapq=$2
			shift
			;;
		-mc|--fb-min-coverage)
			fb_min_cov=$2
			shift
			;;
		-n|--fb-best-alleles)
			nbest=$2
			shift
			;;
		-F|--fb-min-alternate-fraction)
			altfrac=$2
			shift
			;;
		-mvq|--min-variant-quality)
			min_var_qual=$2
			shift
			;;
		-mqa|--min-qual-ao)
			min_qual_ao=$2
			shift
			;;
		-sw|--snp-window)
			snpwindow=$2
			shift
			;;
		-mincov|--minimum-coverage)
			mincov=$2
			shift
			;;
		-annot|--annotate)
			annot=$2
			shift
			;;
		-g|--gff-file)
			gff=$2
			shift
			;;
		-pool|--pooled-sequencing)
			pool=$2
			shift
			;;
		-h|--help)
			echo "$hm"
			exit 1
			;;
		-v|--version)
			echo "0.8.1"
			exit 1
			;;
		*) echo "Unknown parameter passed: $1"
			exit 1
			;;
	esac
	shift
done

cat <<'END_FIGLET'
  _____     ______
 / _ \ \   / / ___|
| | | \ \ / / |  _
| |_| |\ V /| |_| |
 \__\_\ \_/  \____|
             v0.8.1
END_FIGLET

echo "Run started at $(date)"

depends="fastp bwa samtools sambamba freebayes bcftools vcf2fasta vcfstats vcffilter vcftools bedtools bioawk R Rscript" #samtools needs to be at least 1.10!

for i in $depends
do
	if which $i; then
		echo $i found
	else
		echo $i is not found
		echo "Please install $i or specify it in the '$PATH'"
		echo "The pipeline stops now"
		exit 1
	fi
done

if [[ ${#ref_db} -le 1 ]]; then
	echo "Please specify reference genome"
	exit 1
fi

if [[ ${#slist} -le 1 ]]; then
	echo "Please specify samples list"
	exit 1
fi

#check if input files and ref_dbgenom exist

if [[ ! -d ${outdir} ]]; then
	echo "Output directory does not exist"
	exit 1
fi

outdir=`realpath $outdir`

ref_db=`realpath "$ref_db"`

gff=`realpath "$gff"`

inds=$(cut -f 1 "$slist" | sort)

echo "List of samples:"
cat "$slist"

if [[ ${#indir} -le 1 ]]; then
	echo "Please specify input directory of reads to be aligned"
	exit 1
fi

if [[ "$type" == "PE" ]]; then
	bwa index "$ref_db"
	for i in $inds;
	do
		r1=$(find "$indir" -maxdepth 1 -name "${i}*R1*fq.gz" -or -name "${i}*R1*fastq.gz")
		r2=$(find "$indir" -maxdepth 1 -name "${i}*R2*fq.gz" -or -name "${i}*R2*fastq.gz")

		if [[ ! -d "${outdir}/${i}" ]]; then
			mkdir "${outdir}"/"${i}"
		fi

		fsize=$(gzip -l "$r1" | awk 'NR==2 {print $2}')

		if [[ $fsize -gt 0 ]]; then
			echo "##### Filtering ${i} #####"
			fastp -i "$r1" -I "$r2" -o "${outdir}"/"${i}"/"${i}"_R1.fq.gz -O "${outdir}"/"${i}"/"${i}"_R2.fq.gz -5 20 -3 20 -W 10 -M 20 --detect_adapter_for_pe -w "$np" -x -f "$tf1" -F "$tf2" -t "$tt1" -T "$tt2" -l "$ml" -h "${outdir}"/"${i}"/"${i}"_filter.html -j "${outdir}"/"${i}"/"${i}"_filter.json

			echo "##### Aligning ${i} to $ref_db #####"
			bwa mem -t "$np" -k "$bwa_k" -A "$bwa_A" -B "$bwa_B" -O "$bwa_O" -R "@RG\tID:$i\tSM:$i\tPL:Illumina" "$ref_db" "${outdir}"/"${i}"/"${i}"_R1.fq.gz "${outdir}"/"${i}"/"${i}"_R2.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ "$np" |\
			samtools sort -@ "$np" > "${outdir}"/"${i}"/"${i}".bam
		else
			echo "##### Sample file ${i} is empty #####"
			rm -r "${outdir}"/"${i}"
		fi
	done

elif [[ "$type" == "SE" ]]; then
	bwa index "$ref_db"
	for i in $inds;
	do
		r1=$(find "$indir" -maxdepth 1 -name "${i}*R1*fq.gz" -or -name "${i}*R1*fastq.gz")
	
		if [[ ! -d "${outdir}/${i}" ]]; then
			mkdir "${outdir}"/"${i}"
		fi

		fsize=$(gzip -l "$r1" | awk 'NR==2 {print $2}')

		if [[ $fsize -gt 0 ]]; then
			echo "##### Filtering ${i} #####"
			fastp -i "$r1" -o "${outdir}"/"${i}"/"${i}"_R1.fq.gz -5 20 -3 20 -W 10 -M 20 -w "$np" -x -f 10 -t 1 -l 30 -h "${outdir}"/"${i}"/"${i}"_filter.html -j "${outdir}"/"${i}"/"${i}"_filter.json

			echo "##### Aligning ${i} to $ref_db #####"
			bwa mem -t "$np" -k "$bwa_k" -A "$bwa_A" -B "$bwa_B" -O "$bwa_O" -R "@RG\tID:$i\tSM:$i\tPL:Illumina" "$ref_db" "${outdir}"/"${i}"/"${i}"_R1.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ "$np" |\
			samtools sort -@ "$np" > "${outdir}"/"${i}"/"${i}".bam
		else
			echo "##### Sample file ${i} is empty #####"
			rm -r "${outdir}"/"${i}"
		fi
	done

else
	echo "Sequencing type not recognized (SE or PE)"
fi

bioawk -c fastx '{print $name,length($seq)}' "$ref_db" > "${outdir}"/gfile
bedtools makewindows -g "${outdir}"/gfile -w "$clipwindow" -s "$clipstep" > "${outdir}"/bedfile

for i in $inds;
do

	if [[ -f ${outdir}/${i}/${i}.bam ]]; then

		covs=$(samtools coverage "${outdir}"/"${i}"/"${i}".bam | cut -f 7 | tail -1 | awk -v hcm="$hcm" '$2 = $1*hcm')
		highcov=$(echo $covs | cut -f 2 -d " " | cut -f 1 -d ".")
		maxcov=$(samtools depth "${outdir}"/"${i}"/"${i}".bam | cut -f 3 | sort -n | tail -1)

		if [[ $highcov -gt $maxcov ]];then
			echo "##### No clipping of high-coverage regions can be done (threshold is larger than the maximum read depth observed) #####"

			echo "##### Marking duplicate reads of ${i} #####"
			sambamba markdup -p -t "$np" "${outdir}"/"${i}"/"${i}".bam "${outdir}"/"${i}"/"${i}"_markdup.bam

		else
			echo "##### Clipping high coverage regions #####"
			bioawk -c fastx '{print $name,length($seq)}' "$ref_db" > gfile
			samtools index -@ "$np" "${outdir}"/"${i}"/"${i}".bam
			bedtools coverage -a "${outdir}"/bedfile -b "${outdir}"/"${i}"/"${i}".bam > "${outdir}"/"${i}"/"${i}".bedcov
			awk -v md="$highcov" -v OFS="\t" '$4 <= md' "${outdir}"/"${i}"/"${i}".bedcov | bedtools merge > "${outdir}"/"${i}"/"${i}".covbed
			sambamba slice -L "${outdir}"/"${i}"/"${i}".covbed "${outdir}"/"${i}"/"${i}".bam | samtools sort -@ "$np" | samtools view -h -b -@ "$np" > "${outdir}"/"${i}"/"${i}"_clip.bam

			echo "##### Marking duplicate reads of ${i} #####"
			sambamba markdup -p -t "$np" "${outdir}"/"${i}"/"${i}"_clip.bam "${outdir}"/"${i}"/"${i}"_markdup.bam
		fi

		echo "##### Running basic statistics on the alignments of ${i} with samtools #####"	
		
		echo "samtools coverage -m ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_coverage.txt" | parallel -j "$np"
		echo "##### Coverage information for ${i} in histogram (samtools coverage -m) format can be found in ${outdir}/${i}/${i}_coverage.txt #####" 
		echo "samtools coverage ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_coverage.tsv" | parallel -j "$np"
		echo "##### Coverage information for ${i} in tabular format can be found in ${outdir}/${i}/${i}_coverage.tsv #####"
		echo "samtools depth -d 0 ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_depth.tsv" | parallel -j "$np"
		echo "##### Read depth of each position of ${i} in tabular format can be found in ${outdir}/${i}/${i}_depth.tsv #####"
		echo "samtools flagstat ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_flagstat.txt" | parallel -j "$np"
		echo "##### Simple statistics of ${i} (samtools flagstat) can be found in ${outdir}/${i}/${i}_flagstat.txt #####"
		echo "samtools idxstats ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_idxstats.tsv" | parallel -j "$np"
		echo "##### BAM index statistics of ${i} can be found in ${outdir}/${i}/${i}_idxstats.tsv #####"

		if [[ $smooth_coverage == "yes" ]];then

			echo "##### Smoothing coverage with a consecutive window size of $smooth_window and a readcount of $smooth_count within windows #####"
			bedtools makewindows -g "${outdir}"/gfile -w "$smooth_window" > "${outdir}"/cons_bedfile
			mkdir "${outdir}"/"${i}"/loc
			bedtools coverage -sorted -a "${outdir}"/cons_bedfile -b "${outdir}"/"${i}"/"${i}"_markdup.bam |\
				awk -v OFS="\t" -v cov="$smooth_count" '{if ($4 == 0) $8 = 1; else $8 = cov / $4; print $0}' | cut -f 1-3,8 |\
				while read line
				do
					name=$(echo $line | cut -f 1-3 -d " " | sed 's/ /_/g')
					echo "$line" | sed 's/ /\t/g' > "${outdir}"/"${i}"/loc/"${name}"
				done

			cd "${outdir}"/"${i}"/loc/
			beds=$(ls *)
			cd "$cdir"

			mkdir "${outdir}"/"${i}"/temps
			for b in $beds
			do
				prop=$(cut -f 4 "${outdir}"/"${i}"/loc/"${b}")
				sambamba view -L "${outdir}"/"${i}"/loc/"${b}" -h -s "$prop" "${outdir}"/"${i}"/"${i}"_markdup.bam > "${outdir}"/"${i}"/temps/${b}_s.bam 2> /dev/null
			done

			ls "${outdir}"/"${i}"/temps/*s.bam > "${outdir}"/"${i}"/bamlist

			numfiles=$(wc -l "${outdir}"/"${i}"/bamlist | cut -f 1 -d " ")

			limit=$(ulimit -n)

			if [[ $numfiles -ge $limit ]]; then

				newlimit=$(expr "$numfiles+100" | bc)
				ulimit -n "$newlimit"
			fi

			samtools merge -@ "np" -b "${outdir}"/"${i}"/bamlist -f -o "${outdir}"/"${i}"/"${i}"_smooth.bam
			rm "${outdir}"/"${i}"/bamlist
			samtools sort -@ "$np" "${outdir}"/"${i}"/"${i}"_smooth.bam | samtools view -h -b -@ 6 > "${outdir}"/"${i}"/"${i}"_smooth_sort.bam
			mv "${outdir}"/"${i}"/"${i}"_smooth_sort.bam "${outdir}"/"${i}"/"${i}"_smooth.bam
			sambamba markdup -r -t "$np" "${outdir}"/"${i}"/"${i}"_smooth.bam "${outdir}"/"${i}"/"${i}"_markdup.bam 2> /dev/null
			samtools index "${outdir}"/"${i}"/"${i}"_markdup.bam

			rm -r "${outdir}"/"${i}"/temps
			rm -r "${outdir}"/"${i}"/loc
		fi
	fi
done

cd "${outdir}"
cat  <(find ./ -name "*coverage.tsv" | xargs grep -H "^#" | head -n 1 | sed -e 's/\/.*:#/sample_id\t/' -e 's/\.//')  <(find ./ -name "*coverage.tsv" | xargs grep -vH "^#" | sed -e 's/\.\///' -e 's/\/.*:/\t/') > "${outdir}"/coverages.tsv

slist_filt=$(awk -v mincovc="$mincov" '$7 > mincov' "${outdir}"/coverages.tsv | cut -f 1 | grep -v "sample_id")
echo "The following samples cover at least ${mincov} % of the reference genome" #IF ${mincov} IS EMPTY RAISE A WARNING
echo "$slist_filt"
inds=$slist_filt

cd "$cdir"

echo 'args<-commandArgs(TRUE)
x<-read.table(args[1], sep="\t", header=F)
name<-paste(args[1], ".png", sep="")
name<-gsub(".tsv.", ".", name)
png(name, width=3000, height=1000, unit="px")
par(mar=rep(5,4))
plot(x$V3, ylab="Read depth", xlab="Position", type="l", main=args[1], cex.lab=2, cex.axis=2, cex.main=2)
dev.off()

name2<-paste(args[1], "_log.png", sep="")
name2<-gsub(".tsv.", ".", name2)
png(name2, width=3000, height=1000, unit="px")
par(mar=rep(5,4))
plot(log(x$V3), ylab="log(read depth)", xlab="Position", type="l", main=args[1], cex.lab=2, cex.axis=2, cex.main=2)
dev.off()' > "${outdir}"/plot_depth.R

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_depth.tsv ]]; then
		cd "${outdir}"/"${i}"
		Rscript "${outdir}"/plot_depth.R "${i}"_depth.tsv &> /dev/null
		echo "##### Plotting read depth distribution of ${i} #####"
	fi
done

echo 'args<-commandArgs(TRUE)
x<-read.csv(args[1], sep="\t")
name<-paste(args[1], ".pdf", sep="")
name<-gsub(".tsv.", ".", name)
pdf(name, width=20)
plot(x$covbases, xlab=x$sample_id, las=2, pch=19, cex=3, axes=F, main="covered bases (bp)")
axis(2)
box()
plot(x$coverage, xlab=x$sample_id, las=2, pch=19, cex=3, axes=F, main="% of genome covered")
axis(2)
box()
plot(x$numreads, xlab=x$sample_id, las=2, pch=19, cex=3, axes=F, main="no. of reads")
axis(2)
box()
plot(x$meandepth, xlab=x$sample_id, las=2, pch=19, cex=3, axes=F, main="mean read depth (×)")
axis(2)
box()
dev.off()' > "${outdir}"/plot_coverage.R

if [[ -f ${outdir}/coverages.tsv ]]; then
	cd "${outdir}"
	Rscript "${outdir}"/plot_coverage.R coverages.tsv &> /dev/null
	echo "##### Plotting coverage statistics #####"  #COMMENT

fi

echo "##### Summary figure of coverage information for samples can be found in ${outdir}/coverages.pdf #####"

cd "$cdir"

echo "##### Calling variant sites using $np parallel threads with ploidy set to ${ploidy}. Threads are split by sample files. #####"

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
		echo "freebayes -f $ref_db -b ${outdir}/${i}/${i}_markdup.bam -n $nbest --min-coverage $fb_min_cov --ploidy $ploidy -F $altfrac -q $minqual -Q $mismatchqual -m $mapq -w -V -a -j -E -1 > ${outdir}/${i}/${i}_p${ploidy}.vcf"
	fi
done | parallel -j "$np"

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
		vcffilter -f "QUAL > ${min_var_qual} & QUAL / AO > ${min_qual_ao}" "${outdir}"/"${i}"/"${i}"_p"${ploidy}".vcf > "${outdir}"/"${i}"/"${i}"_p"${ploidy}"_filt.vcf
		mv "${outdir}"/"${i}"/"${i}"_p"${ploidy}"_filt.vcf "${outdir}"/"${i}"/"${i}"_p"${ploidy}".vcf
		vcftools --vcf "${outdir}"/"${i}"/"${i}"_p"${ploidy}".vcf  --SNPdensity "${snpwindow}" --out "${outdir}"/"${i}"/"${i}"_SNPdensity 2> /dev/null
	fi
done

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_p${ploidy}.vcf ]]; then
		cd "${outdir}"/"${i}"
		vcfstats "${outdir}"/"${i}"/"${i}"_p"${ploidy}".vcf > "${outdir}"/"${i}"/"${i}"_vcfstats_p"${ploidy}".txt
		echo "##### vcf statistics can be found at ${outdir}/${i}/${i}_vcfstats_p${ploidy}.txt #####"
		vcf2fasta -f "$ref_db" "${i}"_p"${ploidy}".vcf
		bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\t[\t%GT]\n' "${i}"_p"${ploidy}".vcf > "${i}"_p"${ploidy}"_GT.tsv 2> /dev/null
	fi
done

uniqid=$(date | sed -e 's/ /_/g' -e 's/\.//g' -e 's/,//g' | cut -f 1,2,3,5 -d_)

find "${outdir}"/*/*.fa | xargs cat | sed -e 's/\.//' -e 's/://' > "${outdir}/samples_multifasta_${uniqid}.fa"

echo "##### Calling variant sites using $np parallel threads assuming pooled sequencing. Threads are split by sample files. #####"

if [[ $pool == "yes" ]];then

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
			echo "freebayes -f $ref_db -b ${outdir}/${i}/${i}_markdup.bam -n $nbest --min-coverage $fb_min_cov --pooled-continuous -F $altfrac -q $minqual -Q $mismatchqual -m $mapq -w -V -a -j -E -1 > ${outdir}/${i}/${i}_pooled.vcf"
		fi
	done | parallel -j "$np"

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
			vcffilter -f "QUAL > ${min_var_qual} & QUAL / AO > ${min_qual_ao}" "${outdir}"/"${i}"/"${i}"_pooled.vcf > "${outdir}"/"${i}"/"${i}"_pooled_filt.vcf
			mv "${outdir}"/"${i}"/"${i}"_pooled_filt.vcf "${outdir}"/"${i}"/"${i}"_pooled.vcf
		fi
	done

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_pooled.vcf ]]; then
			vcfstats "${outdir}"/"${i}"/"${i}"_pooled.vcf > "${outdir}"/"${i}"/"${i}"_vcfstats_pooled.txt
			echo "##### vcf statistics can be found at ${outdir}/${i}/${i}_vcfstats_pooled.txt #####"
			cd "${outdir}"/"${i}"/
			bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\t[\t%GT]\n' "${i}"_pooled.vcf > "${i}"_pooled_GT.tsv 2> /dev/null
		fi
	done

	echo 'args<-commandArgs(TRUE)
	x<-read.csv(args[1], sep="\t")
	head(x)
	name<-paste(args[1], ".pdf", sep="")
	name<-gsub(".tsv.", ".", name)
	pdf(name, width=20)
	hist(x$X.5.AB, main=args[1], xlab="Allele balance")
	dev.off()' > "${outdir}"/plot_AB.R

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_pooled_GT.tsv ]]; then
			cd "${outdir}"/"${i}"
			cut -f 5 "${i}"_pooled_GT.tsv | perl -pe "s/,/\n/g" > "${i}"_pooled_AB.tsv
			Rscript "${outdir}"/plot_AB.R "${i}"_pooled_AB.tsv &> /dev/null
			echo "##### Plotting AB distribution of ${i} #####"
		fi
	done

	cd "$cdir"

	echo "##### Exporting distribution of sites with more than one allele across samples to ${outdir}/non_haploid_sites.txt #####"
	find "${outdir}" -name *pooled.vcf | xargs grep '0/1:' | cut -f 2 | sort -n | uniq -c > "${outdir}"/non_haploid_sites.txt # instead check for allele balance and export sites with AB > 0

fi

cd "$cdir"

echo "##### Masking of low-depth positions in the consensus fasta files. #####"
for i in $inds
do
	echo "$i"
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]]; then
		bioawk -c fastx '{print $name, length($seq)}' "${outdir}"/"${i}"/"${i}"*.fa > "${outdir}"/"${i}"/"${i}".genomfile
		bwa index "${outdir}"/"${i}"/"${i}"*.fa

		if [[ "$type" == "PE" ]]; then
			echo "##### Realigning ${i} #####"
			samtools sort -n -@ "$np" "${outdir}"/"${i}"/"${i}"_markdup.bam | samtools fastq -@ "$np" -1 "${outdir}"/"${i}"/1.fastq -2 "${outdir}"/"${i}"/2.fastq -s "${outdir}"/"${i}"/s.fastq
			bwa mem -t "$np" -R "@RG\tID:$i\tSM:$i\tPL:Illumina" "${outdir}"/"${i}"/"${i}"*.fa "${outdir}"/"${i}"/1.fastq "${outdir}"/"${i}"/2.fastq 2> /dev/null |\
			samtools view -h -b -u -@ "$np" |\
			samtools sort -@ "$np" > "${outdir}"/"${i}"/realigned.bam
		elif [[ "$type" == "SE" ]]; then
			echo "##### Realigning ${i} #####"
			samtools sort -n -@ "$np" "${outdir}"/"${i}"/"${i}"_markdup.bam | samtools fastq -@ "$np" -0 "${outdir}"/"${i}"/s.fastq
			bwa mem -t "$np" -R "@RG\tID:$i\tSM:$i\tPL:Illumina" "${outdir}"/"${i}"/"${i}"*.fa "${outdir}"/"${i}"/s.fastq 2> /dev/null |\
			samtools view -h -b -u -@ "$np" |\
			samtools sort -@ "$np" > "${outdir}"/"${i}"/realigned.bam
		fi

		echo "##### Finding and masking blocks with < ${fb_min_cov} reads and potential clipped high-coverage regions #####"
		bedtools genomecov -ibam "${outdir}"/"${i}"/realigned.bam -d | awk -v mincov="$fb_min_cov" '$3 < mincov' | awk 'OFS="\t"{print $1, $2-1, $2}' | bedtools merge -i - > "${outdir}"/"${i}"/complement
		bedtools maskfasta -fi "${outdir}"/"${i}"/"${i}"*.fa -bed "${outdir}"/"${i}"/complement -fo "${outdir}"/"${i}"/"${i}"_masked.fasta
	fi
done

find "${outdir}"/*/*masked.fasta | xargs cat | sed -e 's/\.//' -e 's/://' > "${outdir}/samples_multifasta_masked_${uniqid}.fa"

if [ ! -f "$gff" ]; then
	echo 'args<-commandArgs(TRUE)
	x<-read.csv(args[1], sep="\t")
	name<-paste(args[1], ".pdf", sep="")
	name<-gsub(".tsv.", ".", name)
	pdf(name, width=20)
	plot(x$VARIANTS.KB~x$BIN_START, xlab="bin start", ylab="variants/kb", lwd=2, type="l", main=args[1])
	dev.off()' > "${outdir}"/plot_density.R

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_SNPdensity.snpden ]]; then
			cd "${outdir}"/"${i}"
			Rscript "${outdir}"/plot_density.R "${i}"_SNPdensity.snpden &> /dev/null
			echo "##### Plotting SNP density distribution of ${i} #####"  #COMMENT
		fi
	done

	cd "$cdir"

elif [ -f "$gff" ]; then

	echo 'args<-commandArgs(TRUE)
	x<-read.csv(args[1], sep="\t")
	name<-paste(args[1], ".pdf", sep="")
	name<-gsub(".tsv.", ".", name)
	regions<-read.csv(args[2], sep = "\t", header = F)
	pdf(name, width=20)
	plot(x$VARIANTS.KB~x$BIN_START, xlab="bin start", ylab="variants/kb", lwd=2, type="l", ylim=c(0,max(x$VARIANTS.KB+2)), main=args[1])
	rect(xleft=regions$V1, xright = regions$V2, ybottom=max(x$VARIANTS.KB), ytop=max(x$VARIANTS.KB)+1, col=c("#BEBEBE50","#A9A9A950"))
	text(x = (regions$V2+regions$V1)/2, y = max(x$VARIANTS.KB)+1.5, labels = regions$V3, srt=90, cex=0.8)
	dev.off()' > "${outdir}"/plot_density.R

	echo "##### Transfering annotations #####"

	cut -f 3 "$gff" | grep -v "#" | sort | uniq > ${outdir}/feature.types

	for i in $inds
	do
		if [[ -f ${outdir}/${i}/${i}_SNPdensity.snpden ]]; then
			cd "${outdir}"/"${i}"
			liftoff -p "$np" -g "$gff" -mismatch 10000 -a 0.01 -f "${outdir}"/feature.types -o "${i}".gff "${i}"_masked.fasta "$ref_db" &> /dev/null
			grep -v "^#" "${i}".gff | awk '$3 != "CDS"' | awk '$3 != "region"' | cut -f 4,5,9 | sed -e "s/\tID.*Name=/\t/" -e "s/\tID.*gbkey=/\t/" | cut -f 1 -d";" > "${i}"_regions.tsv
			Rscript "${outdir}"/plot_density.R "${i}"_SNPdensity.snpden "${i}"_regions.tsv
			echo "##### Plotting SNP density distribution of ${i} #####"  #COMMENT
			rm -r "${outdir}"/"${i}"/intermediate_files
			rm -r "${outdir}"/"${i}"/*mmi
		fi
	done
	rm "${outdir}"/feature.types
fi

rm "${outdir}"/plot*R
rm "${outdir}"/bedfile
rm "${outdir}"/gfile
find "${outdir}"/*/ -name "complement" | xargs rm
find "${outdir}"/*/ -name "realigned.bam" | xargs rm
find "${outdir}"/*/ -name "*fastq" | xargs rm
find "${outdir}"/*/ -name "*bwt" | xargs rm
find "${outdir}"/*/ -name "*pac" | xargs rm
find "${outdir}"/*/ -name "*amb" | xargs rm
find "${outdir}"/*/ -name "*ann" | xargs rm
find "${outdir}"/*/ -name "*sa" | xargs rm
find "${outdir}"/*/ -name "*.genomfile" | xargs rm

echo "Run ended at $(date)"
#end
