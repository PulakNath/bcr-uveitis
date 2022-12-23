# This script Reads all the cellranger summary metric files for all the samples
# and plots the aggregated data for easy comparison
# Load the below libraries, install them if not already installed
# Author : Vijay Nagarajan PhD, NEI/NIH
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(ggplot2)
library(ggdendro)
library(plotly)
library(htmlwidgets)

##################################### Setup
# Set the input, output paths
directorypath="/data/../TotalSeq/Results/cellranger/"
outputpath="/data/../TotalSeq/Reports/"
sampleinfopath="/data/../TotalSeq/SampleInformation.tab"

# Read the results folders
file_list <- list.files(path=directorypath)
# List result folders
file_list
# Set the column headers
filehead=c("SampleName","Estimated # of Cells","Mean Reads/Cell","Median Genes/Cell","# of Reads","Valid Barcodes","Sequencing Saturation","Q30 Bases in Barcode","Q30 Bases in RNA Read","Q30 Bases in UMI","Reads Mapped to Genome","Reads Mapped Confidently to Genome","Reads Mapped Confid to Intergenic Reg","Reads Mapped Confid to Intronic Reg","Reads Mapped Confid to Exonic Reg","Reads Mapped Confid to Transcriptome","Reads Mapped Antisense to Gene","Fraction Reads in Cells","Total Genes Detected","Median UMI Counts/Cell","Ab: Number of Reads","Ab: Mean Reads/Cell","Ab: Valid Barcodes","Ab: Sequencing Saturation","Ab: Q30 Bases in Barcode","Ab: Q30 Bases in Ab Read","Ab: Q30 Bases in UMI","Ab: Fraction Ab Reads","Ab: Fraction Ab Reads Usable","Ab: Ab Reads Usable/Cell","Ab: Fraction Ab Reads in Aggregate Barcodes","Ab: Fraction Unrecognized Ab","Ab: Ab Reads in Cells","Ab: Median UMIs/Cell (summed over all recognized ab barcodes)")
filehead2=c("Estimated # of Cells","Mean Reads/Cell","Median Genes/Cell","# of Reads","Valid Barcodes","Sequencing Saturation","Q30 Bases in Barcode","Q30 Bases in RNA Read","Q30 Bases in UMI","Reads Mapped to Genome","Reads Mapped Confidently to Genome","Reads Mapped Confid to Intergenic Reg","Reads Mapped Confid to Intronic Reg","Reads Mapped Confid to Exonic Reg","Reads Mapped Confid to Transcriptome","Reads Mapped Antisense to Gene","Fraction Reads in Cells","Total Genes Detected","Median UMI Counts/Cell","Ab: Number of Reads","Ab: Mean Reads/Cell","Ab: Valid Barcodes","Ab: Sequencing Saturation","Ab: Q30 Bases in Barcode","Ab: Q30 Bases in Ab Read","Ab: Q30 Bases in UMI","Ab: Fraction Ab Reads","Ab: Fraction Ab Reads Usable","Ab: Ab Reads Usable/Cell","Ab: Fraction Ab Reads in Aggregate Barcodes","Ab: Fraction Unrecognized Ab","Ab: Ab Reads in Cells","Ab: Median UMIs/Cell (summed over all recognized ab barcodes)")
# Initiate an empty dataframe
filealldata = data.frame()

##################################### Extract and map short names
# Read sample information file
sampleinformation=read.csv(sampleinfopath,sep="\t")
# Extract shortname columns
samplenametype=unique(paste(sampleinformation$submitted_name,sampleinformation$ShortNameType,sep=","))
# Convert to dataframe
samplenametype=as.data.frame(samplenametype)
# Initiate emply lists to save hc and uv sample names
samplehclist=list()
sampleuvlist=list()

##################################### Extract and aggregate summary metrics
# Read individual summary metric file
for(i in 1:length(file_list))
{
  # File path
  filename=paste(directorypath,file_list[i],"/cellranger_output/metrics_summary.csv",sep = "")
  # Store sample name
  filesample=unlist(file_list[i])
  # Identify and store short sample name
  filesampleoneline=samplenametype[grep(filesample,samplenametype[,1]), ]
  filesampleoneline=unlist(as.vector(strsplit(filesampleoneline, ',')))
  filesampleoneline=filesampleoneline[2]
  print(filesampleoneline)
  filesample=filesampleoneline
  # Create HC name list
  if(grepl("HC",filesample))
  {
    samplehclist<-append(samplehclist,filesample)
  }
  # Create UV name list
  if(grepl("UV",filesample))
  {
    sampleuvlist<-append(sampleuvlist,filesample)
  }
  # Read summary metric data
  filedata <- read_csv(filename, col_types = cols(`Valid Barcodes` = col_number(), `Sequencing Saturation` = col_number(), `Q30 Bases in Barcode` = col_number(), `Q30 Bases in RNA Read` = col_number(), `Q30 Bases in UMI` = col_number(), `Reads Mapped to Genome` = col_number(), `Reads Mapped Confidently to Genome` = col_number(), `Reads Mapped Confidently to Intergenic Regions` = col_number(), `Reads Mapped Confidently to Intronic Regions` = col_number(), `Reads Mapped Confidently to Exonic Regions` = col_number(), `Reads Mapped Confidently to Transcriptome` = col_number(), `Reads Mapped Antisense to Gene` = col_number(), `Fraction Reads in Cells` = col_number(), `Antibody: Valid Barcodes` = col_number(), `Antibody: Sequencing Saturation` = col_number(), `Antibody: Q30 Bases in Barcode` = col_number(), `Antibody: Q30 Bases in Antibody Read` = col_number(), `Antibody: Q30 Bases in UMI` = col_number(), `Antibody: Fraction Antibody Reads` = col_number(), `Antibody: Fraction Antibody Reads Usable` = col_number(), `Antibody: Fraction Antibody Reads in Aggregate Barcodes` = col_number(), `Antibody: Fraction Unrecognized Antibody` = col_number(), `Antibody: Antibody Reads in Cells` = col_number()))
  # Aggregate summary metric data in to filealldata variable
  if(length(filedata[2,])>=33)
  {
    filedataframe=as.data.frame(filedata[1,])
    filedatanumbers=as.numeric(filedataframe[1,])
    fileonerow=c(filesample,filedatanumbers)
    filealldata=rbind(filealldata,c(filesample,filedatanumbers))
  }
}

##################################### Plot and save interactive multipanel image as html file
# Create list for ordering hc and uv samples
samplehclist=unlist(samplehclist)
sampleuvlist=unlist(sampleuvlist)
sampletypelist=append(samplehclist,sampleuvlist)
print(sampletypelist)

# Add column headers
filealldata
forplot=filealldata
colnames(forplot)=filehead

# Exclude samplename column and include all other columns for plotting
dfm <- pivot_longer(forplot, -SampleName, names_to="variable", values_to="value")
intplot=ggplot(dfm, aes(x = factor(SampleName, levels=sampletypelist),y = as.numeric(value))) +
  geom_bar(aes(fill = variable), position="stack", stat="identity")+
  facet_wrap(~variable, scales="free_y")+
  theme(legend.position = "none",
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 90),
    strip.text = element_text(size=8),
    axis.title.x = element_blank())
# Generate the plot and save it as html in the "Reports" folder
l <- plotly::ggplotly(intplot)
htmlwidgets::saveWidget(as_widget(l), paste(outputpath,"cellrangerReport.html",sep=""))

# Cluster sample metrics and generate a dendrogram
nohead=filealldata[,-1]
temp2=filealldata[,1]
forheatmap <- sapply( nohead, as.numeric )
forheatmap
# Add column headers
rownames(forheatmap)=temp2
colnames(forheatmap)=filehead2
# Convert data to matrix
data=as.matrix(forheatmap)
# Calculate the vector similarity using euclidean distance method
dd <- dist(scale(data), method = "euclidean")
# Draw the dendrogram
hc <- hclust(dd, method = "ward.D2")
# Save the dendrogram as an image in the "Reports" folder
png(paste(outputpath,"cellrangerClusterReport.png",sep=""),res=300)
ggdendrogram(hc)
dev.off()
