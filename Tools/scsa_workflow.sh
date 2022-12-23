# The code in this script couild be used to generate scsa annotations
# Author : Vijay Nagarajan PhD, NEI/NIH
# Load python module
module load python/3.8

# Downgrade pandas
pip install --user pandas==1.2.4

# Change to scsa folder
cd /data/../TotalSeq/Tools/SCSA
cp /data/../TotalSeq/Results/markers/markersall.txt /data/../TotalSeq/Results/markers/markersall_forscsa.txt
# Edit markersall_forscsa to change column name 'avg_log2FC' to avg_logFC

# Run scsa with blood reference
python SCSA.py -d whole.db -s seurat -i /data/../TotalSeq/Results/markers/markersall_forscsa.txt -k Blood -g Human -E -p 0.05 -f 2 -m txt -o /data/../TotalSeq/Results/annotation/markersallscsa_dec2022_blood.txt

# Run scsa with all reference
python SCSA.py -d whole.db -s seurat -i /data/../TotalSeq/Results/markers/markersall_forscsa.txt -k All -g Human -E -p 0.05 -f 2 -m txt -o /data/../TotalSeq/Results/annotation/markersallscsa_dec2022_all.txt
