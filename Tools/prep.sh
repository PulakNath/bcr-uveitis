#!/bin/bash
# Author : Vijay Nagarajan PhD, NEI/NIH
# This scripts takes a sampleinformation tab delimited files and samples name, creates samples folder with library.csv and swarm file

# Sample Names from SampleInformation.tab file
samples="NS3R189BTS NS7R65BBTS"

# Swarm file
echo "# This swarm file runs cellranger for all totalseq samples" > totalseqall.swarm

# Output results directory
outputdirectory="/data/../TotalSeq/Results/cellranger/"

# RawData directory
inputdirectory="/data/../TotalSeq/RawData/"

# Sampleinfo path
sampleinfofile="/data/../TotalSeq/SampleInformation.tab"
TAB=$'\t'

# Read one sample at a time and extract relevant info from sampleinfo file
for s in ${samples};
	do
	#echo "${s}"
	mkdir "${outputdirectory}${s}"
	filedirectory=`grep  "${s}" "${sampleinfofile}" | cut -f3 | uniq`
	samplenames=`grep  "${s}" "${sampleinfofile}" | cut -f1 | uniq`
	#echo "${filedirectory}"
	#echo "${samplenames}"
	echo -e "fastqs,sample,library_type" > "${outputdirectory}${s}/library.csv"

## Prepare library files and move it in to the corresponding sample folders
	for samplename in ${samplenames};
		do
		sampletype=`grep  "${samplename}" "${sampleinfofile}" | cut -f12 | uniq`
		#echo "${samplename}"
		#echo "${sampletype}"
		echo "${outputdirectory}${s}/,${samplename},${sampletype}" >> "${outputdirectory}${s}/library.csv"
		done

## Create sample directories in results folder and link appropriate raw data files to corresponding sample folders
	#echo "${outputdirectory}${s}/"
	grep  "${s}" "${sampleinfofile}" | cut -f2 > tempfastqfiles
	for fastqfile in $(cat tempfastqfiles);
		do
		#echo "${fastqfile}"
		ln -s "${filedirectory}${fastqfile}" "${outputdirectory}${s}/${fastqfile}"
		done

## Prepare swarm file
	echo "cellranger count --id=${s} --libraries=${outputdirectory}${s}/library.csv --transcriptome=/fdb/cellranger/refdata-gex-GRCh38-2020-A --feature-ref=/data/../TotalSeq/Tools/feature_reference.csv --expect-cells=10000 --localcores=28 --localmem=64 --chemistry=SC3Pv3" >> totalseqall.swarm

## Printout final contents for visual check
	ls -lh "${outputdirectory}${s}"
	cat "${outputdirectory}${s}/library.csv"
	grep  "${s}" "${sampleinfofile}" | cut -f1,2,3,4,12

## Change folder permission for group read/write access
	chmod -R 770 "${outputdirectory}${s}"
	done

cat totalseqall.swarm
rm tempfastqfiles
