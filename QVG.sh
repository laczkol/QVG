#!/bin/bash
set -euo errexit -euo pipefail -euo nounset

echo "Run started at $(date)"

cdir=$(pwd)
indir=$(pwd)
outdir=$(pwd)
ref_db=" "

np=1
map="yes"
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
		-r|--ref_dberence-genome) #ref_dberence genome in fasta, MANDATORY
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
		-trim_front1|--trim_front1)
			tf1=$2
			shift
			;;
		-trim_front2|--trim_front2)
			tf2=$2
			shift
			;;
		-trim_tail1|--trim_tail1)
			tt1=$2
			shift
			;;
		-trim_tail2|--trim_tail2)
			tt2=$2
			shift
			;;
		-minlen|--min-read-length)
			ml=$2
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
		*) echo "Unknown parameter passed: $1"
			exit 1
			;;
	esac
	shift
done



depends="fastp bwa samtools sambamba freebayes bcftools vcf2fasta vcfstats bedtools" #samtools needs to be at least 1.10!

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
	echo "Please specify ref_dberence genome"
	exit 1
fi

if [[ ${#slist} -le 1 ]]; then
	echo "Please specify slist"
	exit 1
fi

#check if input files and ref_dbgenom exist

if [[ ! -d ${outdir} ]]; then
	echo "Output directory does not exist"
	exit 1
fi

outdir=`realpath $outdir`

ref_db=`realpath $ref_db`

inds=`cut -f 1 $slist | sort`

echo "List of samples:"
cat $slist

if [[ ${#indir} -le 1 ]]; then
	echo "Please specify input directory of reads to be aligned"
	exit 1
fi

if [[ "$type" == "PE" ]]; then
	bwa index $ref_db
	for i in $inds;
	do
		r1=`find $indir -maxdepth 1 -name "${i}*R1*fq.gz" -or -name "${i}*R1*fastq.gz"`
		r2=`find $indir -maxdepth 1 -name "${i}*R2*fq.gz" -or -name "${i}*R2*fastq.gz"`

		if [[ ! -d "${outdir}/${i}" ]]; then
			mkdir ${outdir}/${i}
		fi

		fsize=`gzip -l $r1 | awk 'NR==2 {print $2}'`

		if [[ $fsize -gt 0 ]]; then
			echo "##### Filtering ${i} #####"
			fastp -i $r1 -I $r2 -o ${outdir}/${i}/${i}_R1.fq.gz -O ${outdir}/${i}/${i}_R2.fq.gz -5 20 -3 20 -W 10 -M 20 --detect_adapter_for_pe -w $np -x -f $tf1 -F $tf2 -t $tt1 -T $tt2 -l $ml -h ${outdir}/${i}/${i}_filter.html -j ${outdir}/${i}/${i}_filter.json

			echo "##### Aligning ${i} to $ref_db #####"
			bwa mem -t $np -R "@RG\tID:$i\tSM:$i\tPL:Illumina" $ref_db ${outdir}/${i}/${i}_R1.fq.gz ${outdir}/${i}/${i}_R2.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ $np |\
			samtools sort -@ $np > ${outdir}/${i}/${i}.bam
		else
			echo "##### Sample file ${i} is empty #####"
			rm -r ${outdir}/${i}
		fi
	done

elif [[ "$type" == "SE" ]]; then
	bwa index $ref_db
	for i in $inds;
	do
		r1=`find $indir -maxdepth 1 -name "${i}*R1*fq.gz" -or -name "${i}*R1*fastq.gz"`
	
		if [[ ! -d "${outdir}/${i}" ]]; then
			mkdir ${outdir}/${i}
		fi

		fsize=`gzip -l $r1 | awk 'NR==2 {print $2}'`

		if [[ $fsize -gt 0 ]]; then
			echo "##### Filtering ${i} #####"
			fastp -i $r1 -o ${outdir}/${i}/${i}_R1.fq.gz -5 20 -3 20 -W 10 -M 20 -w $np -x -f 10 -t 1 -l 30 -h ${outdir}/${i}/${i}_filter.html -j ${outdir}/${i}/${i}_filter.json

			echo "##### Aligning ${i} to $ref_db #####"
			bwa mem -t $np -R "@RG\tID:$i\tSM:$i\tPL:Illumina" $ref_db ${outdir}/${i}/${i}_R1.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ $np |\
			samtools sort -@ $np > ${outdir}/${i}/${i}.bam
		else
			echo "##### Sample file ${i} is empty #####"
			rm -r ${outdir}/${i}
		fi
	done

else
	echo "Sequencing type not recognized (SE or PE)"
fi

for i in $inds;
do

	if [[ -f ${outdir}/${i}/${i}.bam ]]; then
		echo "##### Marking duplicate reads of ${i} #####"
		sambamba markdup -p -t $np ${outdir}/${i}/${i}.bam ${outdir}/${i}/${i}_markdup.bam

		echo "##### Running basic statistics on the alignments of ${i} with samtools #####"	
		
		echo "samtools coverage -m ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_coverage.txt" | parallel -j $np
		echo "##### Coverage information for ${i} in histogram (samtools coverage -m) format can be found in ${outdir}/${i}/${i}_coverage.txt #####" 
		echo "samtools coverage ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_coverage.tsv" | parallel -j $np
		echo "##### Coverage information for ${i} in tabular format can be found in ${outdir}/${i}/${i}_coverage.tsv #####"
		echo "samtools depth -d 0 ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_depth.tsv" | parallel -j $np
		echo "##### Read depth of each position of ${i} in tabular format can be found in ${outdir}/${i}/${i}_depth.tsv #####"
		echo "samtools flagstat ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_flagstat.txt" | parallel -j $np
		echo "##### Simple statistics of ${i} (samtools flagstat) can be found in ${outdir}/${i}/${i}_flagstat.txt #####"
		echo "samtools idxstats ${outdir}/${i}/${i}_markdup.bam > ${outdir}/${i}/${i}_idxstats.tsv" | parallel -j $np
		echo "##### BAM index statistics of ${i} can be found in ${outdir}/${i}/${i}_idxstats.tsv #####"
	fi
done


echo "##### Calling variant sites using $np parallel threads with ploidy set to ${ploidy}. Threads are split by sample files. #####"

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
		echo "freebayes -f $ref_db -b ${outdir}/${i}/${i}_markdup.bam -n $nbest --min-coverage $fb_min_cov --ploidy $ploidy -F $altfrac -q $minqual -Q $mismatchqual -m $mapq -w -V -a -j -E -1 > ${outdir}/${i}/${i}_p${ploidy}.vcf"
	fi
done | parallel -j $np

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_p${ploidy}.vcf ]]; then
		vcfstats ${outdir}/${i}/${i}_p${ploidy}.vcf > ${outdir}/${i}/${i}_vcfstats_p${ploidy}.txt
		echo "##### vcf statistics can be found at ${outdir}/${i}/${i}_vcfstats_p${ploidy}.txt #####"
		cd ${outdir}/${i}/
		vcf2fasta -f $ref_db ${i}_p${ploidy}.vcf
		bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\t[\t%GT]\n' ${i}_p${ploidy}.vcf > ${i}_p${ploidy}_GT.tsv 2> /dev/null
	fi
done

cd $cdir

uniqid=`date | sed -e 's/ /_/g' -e 's/\.//g' -e 's/,//g' | cut -f 1,2,3,5 -d_`

find ${outdir}/*/*.fa | xargs cat | sed -e 's/\.//' -e 's/://' > "${outdir}/samples_multifasta_${uniqid}.fa"

echo "##### Calling variant sites using $np parallel threads assuming pooled sequencing. Threads are split by sample files. #####"

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]];then
		echo "freebayes -f $ref_db -b ${outdir}/${i}/${i}_markdup.bam -n $nbest --min-coverage $fb_min_cov -K -F $altfrac -q $minqual -Q $mismatchqual -m $mapq -w -V -a -j -E -1 > ${outdir}/${i}/${i}_pooled.vcf"
	fi	
done | parallel -j $np

for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_pooled.vcf ]]; then
		vcfstats ${outdir}/${i}/${i}_pooled.vcf > ${outdir}/${i}/${i}_vcfstats_pooled.txt
		echo "##### vcf statistics can be found at ${outdir}/${i}/${i}_vcfstats_pooled.txt #####"
		cd ${outdir}/${i}/
		bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\t[\t%GT]\n' ${i}_pooled.vcf > ${i}_pooled_GT.tsv 2> /dev/null
	fi
done

cd $cdir


#annotate vcf with gff3 values


#subset fasta genomes by bed region (that is the sequenced region of the genome)

echo "##### Masking of unsequenced positions #####"
for i in $inds
do
	if [[ -f ${outdir}/${i}/${i}_markdup.bam ]]; then
		samtools faidx ${outdir}/${i}/${i}*.fa #should be the only fasta; if not unexpected behavior may follow, watch out!
		cut -f 1-2 ${outdir}/${i}/${i}*.fa.fai > ${outdir}/${i}/${i}.genomfile
		bwa index ${outdir}/${i}/${i}*.fa

		if [[ "$type" == "PE" ]]; then
			echo "##### Realigning ${i} #####"
			bwa mem -t $np -R "@RG\tID:$i\tSM:$i\tPL:Illumina" ${outdir}/${i}/${i}*.fa ${outdir}/${i}/${i}_R1.fq.gz ${outdir}/${i}/${i}_R2.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ $np |\
			samtools sort -@ $np > ${outdir}/${i}/realigned.bam
		elif [[ "$type" == "SE" ]]; then
			echo "##### Realigning ${i} #####"
			bwa mem -t $np -R "@RG\tID:$i\tSM:$i\tPL:Illumina" ${outdir}/${i}/${i}*.fa ${outdir}/${i}/${i}_R1.fq.gz 2> /dev/null |\
			samtools view -h -b -u -@ $np |\
			samtools sort -@ $np > ${outdir}/${i}/realigned.bam
		fi

		echo "##### Finding blocks with zero reads #####"
		bedtools bamtobed -i ${outdir}/${i}/realigned.bam | bedtools merge -i - > ${outdir}/${i}/realigned.bed
		bedtools complement -i ${outdir}/${i}/realigned.bed -g ${outdir}/${i}/${i}.genomfile > ${outdir}/${i}/complement
		bedtools maskfasta -fi ${outdir}/${i}/${i}*.fa -bed ${outdir}/${i}/complement -fo ${outdir}/${i}/${i}_masked.fasta
	fi
done

find ${outdir}/*/*masked.fasta | xargs cat | sed -e 's/\.//' -e 's/://' > "${outdir}/samples_multifasta_masked_${uniqid}.fa"

echo "Run ended at $(date)"
#end
