#!/usr/bin/bash

# set -eo pipefail

prefix="all_samples"
workdir="/data/local/agosdsc/projects/mrg1_nanopore/flair_analyzed/pooled_hom_mrg1_v3"
manifest_file="/home/agosdsc/projects/mrg1_nanopore/as_detection_flair/reads_manifests_pooled_hom_mrg1_v3.tsv"
fastq_dir="/data/local/agosdsc/projects/mrg1_nanopore/results/pooled_hom_MRG1/merged_fastq/"
genome_file="/data/local/agosdsc/projects/mrg1_nanopore/refGenome/WBcel235/Caenorhabditis_elegans.WBcel235.dna.fa"

gtf_file="/data/local/agosdsc/projects/mrg1_nanopore/refGenome/WBcel235/Caenorhabditis_elegans.WBcel235.97.gtf"
short_reads_dir="/fast/AG_Akalin/agosdsc/projects/GASSER_mrg1_rnaseq/pigx_rnaseq_results/mapped_reads/"
short_reads_files=$(find ${short_reads_dir}*.sortedByCoord.out.bam | head -n 2 | tr '\n' ',' | sed 's/.$//' )
threads=32
minSRsupport=1


subset_control="N2"
subset_case="mrg1"
skip_samples="mrg11b"

fastq_files_space=$(cat ${manifest_file} | cut -f4 | tr '\n' ' ' )
fastq_files_comma=$(cat ${manifest_file} | cut -f4 | tr '\n' ',' | sed 's/.$//' )

flair_src="./lib/flair"

if [ -n $prefix ]; then
	prefix="${prefix}_"
fi

if [ ! -d ${workdir} ]; then 
	mkdir ${workdir} 
fi

## initial mapping
echo "echo initial mapping"
dir="${workdir}/align"
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/${prefix}flair.aligned.bam" ]; then
	cmd="python ${flair_src}/flair.py align -g ${genome_file} \
		-r ${fastq_files_space} -t ${threads} \
		-o ${dir}/${prefix}flair.aligned"
	echo "$cmd" 
else echo "echo skipping"
fi

## junction refinement 
echo "echo junction refinement" 
dir=${workdir}/correct
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/shortread_junctions_min${minSRsupport}.bed" ];then
	cmd="python ${flair_src}/bin/junctions_from_sam.py \
		-s ${short_reads_files} \
		-n ${dir}/shortread; \
		awk '"'$5 >= '${minSRsupport}' {print $0}'"' \
		${dir}/shortread_junctions.bed > \
		${dir}/shortread_junctions_min${minSRsupport}.bed
		"
	echo "$cmd"
else echo "echo skipping"
fi

## correct mapping

if [ ! -f "${genome_file}.chromsizes" ]; then
	cmd="samtools faidx ${genome_file}; cut -f1,2 ${genome_file}.fai > ${genome_file}.chromsizes"
	echo "$cmd"
fi

echo "echo correct mapping"
dir=${workdir}/correct_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/${prefix}flair_all_corrected.bed" ]; then
	cmd="python ${flair_src}/flair.py correct -c ${genome_file}.chromsizes	\
		-q ${workdir}/align/${prefix}flair.aligned.bed \
		-g ${genome_file} -f ${gtf_file} \
		-t ${threads} -o ${dir}/${prefix}flair"
	echo "$cmd"
else echo "echo skipping"
fi

## collapse isoforms
echo "echo collapse isoforms"
dir=${workdir}/collapse_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/${prefix}flair.isoforms.fa" ]; then
	cmd="python ${flair_src}/flair.py collapse -g ${genome_file} \
	   -r ${fastq_files_comma} -t ${threads}	\
	   -q ${workdir}/correct/${prefix}flair_all_corrected.psl \
	   -f ${gtf_file} -o ${dir}/${prefix}flair"
	echo "$cmd"
else echo "echo skipping"
fi

## quantify new annotation
echo "echo quantify new annotation"
dir=${workdir}/quantify_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/count_matrix.tsv" ]; then
cmd="python ${flair_src}/flair.py quantify \
	-r ${manifest_file} -t ${threads}	\
	-i ${workdir}/collapse/${prefix}flair.isoforms.fa \
	--tpm -o ${dir}/count_matrix.tsv"
echo "$cmd"
else echo "echo skipping"
fi

##subset count matrix
echo "echo subsetting count matrix"
dir=${workdir}/quantify_noShortReads
subset_prefix="${subset_case}vs${subset_control}"
if [ ! -f "${dir}/count_matrix.${subset_prefix}.tsv" ] && [ -f "${dir}/count_matrix.tsv" ]; then
	subset_skip_cols="$(cut -f1 $manifest_file | grep -e $subset_control  -e $subset_case -n -v | cut -f1 -d ':' )"
	subset_skip_cols="$subset_skip_cols $(cut -f1 $manifest_file | grep -e $skip_samples -n | cut -f1 -d ':' )"
	subset_skip_cols="$( for item in ${subset_skip_cols[*]}; do expr $item + 1; echo; done )"
	subset_skip_cols="$(echo $subset_skip_cols | sed  's/ /,/g')"
	cmd="cut -f$subset_skip_cols --complement ${dir}/count_matrix.tsv > ${dir}/count_matrix.${subset_prefix}.tsv"
	echo "$cmd"
else echo "echo skipping"
fi


# ## differential analysis
echo "echo differential analysis"
dir=${workdir}/diffExp_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -d "${dir}/${subset_prefix}" ] && [ ! -f "${dir}/count_matrix.${subset_prefix}.tsv" ]; then
	cmd="python ${flair_src}/flair.py diffExp \
		-q ${workdir}/quantify_noShortReads/count_matrix.${subset_prefix}.tsv \
		-o ${dir}/${subset_prefix}	-t ${threads}"
	echo "$cmd"
else echo "echo skipping"
fi

# ## differential splicing
echo "echo differential splicing"
dir=${workdir}/diffSplice_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/${subset_prefix}/${subset_prefix}.stderr.txt" ]; then
	cmd="mkdir ${dir}/${subset_prefix}; \
		python ${flair_src}/flair.py diffSplice \
		-i ${workdir}/collapse/${prefix}flair.isoforms.psl \
		-q ${workdir}/quantify_noShortReads/count_matrix.${subset_prefix}.tsv \
		--test -t ${threads} \
		--drim1 3 \
		--drim3 10 \
		-o ${dir}/${subset_prefix}/${subset_prefix}" 
	echo "$cmd"
else echo "echo skipping"
fi

# productivity 
echo "echo predict productivity"
dir=${workdir}/productivity_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/productivity.bed" ]; then
	cmd="python ${flair_src}/bin/predictProductivity.py \
		-i ${workdir}/collapse/${prefix}flair.isoforms.psl \
		-g ${gtf_file} \
		-f ${genome_file} \
		--longestORF \
		> ${dir}/productivity.bed"
	echo "$cmd"
else echo "echo skipping"
fi

# intron retention
echo "echo mark intron retention"
dir=${workdir}/intronRetention_noShortReads
if [ ! -d ${dir} ]; then 
	mkdir ${dir} 
fi
if [ ! -f "${dir}/${prefix}out_coords.txt" ]; then
	cmd="python ${flair_src}/bin/mark_intron_retention.py \
		${workdir}/collapse/${prefix}flair.isoforms.psl \
		${dir}/${prefix}out_isoforms.psl \
		${dir}/${prefix}out_coords.txt"
	echo "$cmd"
else echo "echo skipping"
fi

