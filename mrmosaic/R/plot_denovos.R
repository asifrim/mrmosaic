#!/usr/bin/env /software/bin/Rscript
library(ggplot2)
library(argparse)

VERSION=0.1

combine_data <- function(child_file, dad_file, mum_file){
  dtmp <- read.delim(child_file,header=T)
  dtmp$sample <- "child"
  d <- dtmp
  dtmp <- read.delim(dad_file,header=T)
  dtmp$sample <- "father"
  d <- rbind(d,dtmp)
  dtmp <- read.delim(mum_file,header=T)
  dtmp$sample <- "mother"
  d <- rbind(d,dtmp)
  return(d)
}

plot_trio_chromosome <- function(data, denovos, chromosome,outputdir,sample_id){
  calls_chrom <- denovos[denovos$Chrom == chromosome,]
  
  ds <- data[data$Chr == chromosome,]
  ds$panel <- "Log.R.Ratio"
  ds[ds$Log.R.Ratio > 1,"Log.R.Ratio"] <- 1 
  ds[ds$Log.R.Ratio < -1,"Log.R.Ratio"] <- -1 
  dshet <- ds[ds$GType == "AB",]
  dshet$panel <- "B-Allele Freq"
  
  calls_chrom_adm <- calls_chrom
  calls_chrom_baf <- calls_chrom
  calls_chrom_adm$panel <- "Log.R.Ratio"
  calls_chrom_baf$panel <- "B-Allele Freq"
  
  
  q <- ggplot(ds) + facet_grid(panel~sample,scale="free") + theme_bw() + theme(strip.background= element_rect(fill = "white",color = "white"))
  q <- q + geom_rect(data=calls_chrom_adm, mapping=aes(xmin=IniPos, xmax=EndPos, ymin = -1, ymax=1, fill=type) )
  q <- q + geom_rect(data=calls_chrom_baf, mapping=aes(xmin=IniPos, xmax=EndPos, ymin = 0, ymax=1, fill=type) )
  q <- q + geom_point(data=ds, aes(x=Position, y=Log.R.Ratio),size=0.5)
  q <- q + geom_point(data=dshet, aes(x=Position, y=B.Allele.Freq), color="red",size=0.5) 
  q <- q + ylab("")
  
  output_path = paste(outputdir,"/",sample_id,"_","chr",chromosome,".pdf",sep="")
  ggsave(filename=output_path,plot=q, width=16, height=10)
  
  
}

plot_trio_chromosome_adm <- function(data, denovos, chromosome, outputdir, sample_id){
  calls_chrom <- denovos[denovos$Chrom == chromosome,]
  
  ds <- data[data$Chr == chromosome,]
  ds$panel <- "ADM3"
  ds[ds$ADM3 > 10,"ADM3"] <- 10 
  ds[ds$ADM3 < -10,"ADM3"] <- -10 
  dshet <- ds[ds$GType == "AB",]
  dshet$panel <- "B-Allele Freq"
  
  calls_chrom_adm <- calls_chrom
  calls_chrom_baf <- calls_chrom
  calls_chrom_adm$panel <- "ADM3"
  calls_chrom_baf$panel <- "B-Allele Freq"
  
  
  q <- ggplot(ds) + facet_grid(panel~sample,scale="free") + theme_bw() + theme(strip.background= element_rect(fill = "white",color = "white"))
  q <- q + geom_rect(data=calls_chrom_adm, mapping=aes(xmin=IniPos, xmax=EndPos, ymin = -10, ymax=10, fill=type) )
  q <- q + geom_rect(data=calls_chrom_baf, mapping=aes(xmin=IniPos, xmax=EndPos, ymin = 0, ymax=1, fill=type) )
  q <- q + geom_point(data=ds, aes(x=Position, y=ADM3),size=0.5)
  q <- q + geom_point(data=dshet, aes(x=Position, y=B.Allele.Freq), color="red",size=0.5) 
  q <- q + ylab("")
  
  output_path = paste(outputdir,"/",sample_id,"_","chr",chromosome,".pdf",sep="")
  ggsave(filename=output_path,plot=q, width=16, height=10)
  
}

plot_denovos <- function(child_file, dad_file, mum_file, denovo_file, output_dir, sample_id){
  d <- combine_data(child_file, dad_file, mum_file)
  denovos <- read.table(denovo_file,header=T)
  denovos <- denovos[denovos$inheritance == "denovo",]
  denovos$size <- denovos$EndPos - denovos$IniPos
  denovos <- denovos[denovos$size > 5000000,]
  chromosomes <- unique(denovos$Chrom)
  
  for(chrom in chromosomes){
    plot_trio_chromosome_adm(d,denovos,chrom,output_dir,sample_id)
  }
  
}


parser <- argparse::ArgumentParser(description=paste("MrMosaic Trio De Novo Plotter (version ",VERSION,"):"))
parser$add_argument("--child", help = "Path to child *.adm3.txt file")
parser$add_argument("--dad", help = "Path to dad *.adm3.txt file")
parser$add_argument("--mum", help = "Path to mum *.adm3.txt file")
parser$add_argument("--denovo", help = "Output file of the trio caller script")
parser$add_argument("--prefix", help = "Prefix for output files")
parser$add_argument("--outputdir", help = "Path to directory to output to")

args <- parser$parse_args()

child_file <- args$child
dad_file <- args$dad
mum_file <- args$mum
denovo_file <- args$denovo

output_dir <- args$outputdir

sample_id <- args$prefix


plot_denovos(child_file, dad_file, mum_file, denovo_file, output_dir, sample_id)
