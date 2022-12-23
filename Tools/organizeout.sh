#!/bin/bash
# This scripts is run after cellranger run, to move output files to the results folder
# Author : Vijay Nagarajan PhD, NEI/NIH

# Sample Names
samples="NS3R189BTS NS7R65BBTS"

# Output results directory
outputdirectory="/data/../TotalSeq/Results/cellranger/"

# Input results directory
inputdirectory="/data/../TotalSeq/Tools/"

TAB=$'\t'

# move swarm logs
mkdir -p /data/../TotalSeq/Tools/swarmLogs
mv /data/../TotalSeq/Tools/swarm_* swarmLogs/.

# Read one sample at a time and extract relevant info from sampleinfo file
for s in ${samples};
	do
	mkdir -p "${outputdirectory}${s}/cellranger_output"
	mkdir -p "${outputdirectory}${s}/cellranger_output/cellranger_log"
	remove="../Results/cellranger/${s}/*.fastq.gz"
	rm ${remove}
	echo ${s}
	#rm -rf "${inputdirectory}${s}/outs"
	#cp -R "${inputdirectory}${s}/*" "${outputdirectory}${s}/cellranger_output/cellranger_log/."
	in="/data/../TotalSeq/Tools/${s}/outs/*"
	out="../Results/cellranger/${s}/cellranger_output/."
	in2="/data/../TotalSeq/Tools/${s}/*"
	out2="../Results/cellranger/${s}/cellranger_output/cellranger_log/."

	# move outs folder
	mv ${in} ${out}
	rm -rf "/data/../TotalSeq/Tools/${s}/outs"
	mv ${in2} ${out2}
	rm -rf "/data/../TotalSeq/Tools/${s}"

	## Change folder permission for group read/write access
	chmod -R 770 "${outputdirectory}${s}"

	done
