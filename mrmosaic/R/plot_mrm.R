#!/usr/bin/env Rscript
VERSION=0.1

library(argparse)
library(yaml)
library(gada)
library(zoo)
library(parallel)
library(mrmosaic)


####################################################################################################
## ARGUMENT PARSER


#Parser
parser <- argparse::ArgumentParser(description=paste("MrMosaic plotter (version ",VERSION,"):"))
parser$add_argument("--child_mrm", help = "Path to child mrm file (output of cnv_baf_data)")
parser$add_argument("--dad_mrm", help = "Path to dad mrm file (output of cnv_baf_data)")
parser$add_argument("--mum_mrm", help = "Path to mum mrm file (output of cnv_baf_data)")
parser$add_argument("--chrom", help = "Comma-separated list of chromosomes to plot")
parser$add_argument("--calls", help = "Path to file containing MrMosaic calls")
parser$add_argument("--config",help= "Path to config file")
args <- parser$parse_args()

if (is.null(args$config) || is.null(args$input)){
	parser$print_help()
	stop("Error: You must provide an --input file and a --config file")
}

config <- yaml::yaml.load_file("/nfs/users/nfs_a/as33/Projects/MrMosaic/config.yaml")
rd_data_dir = config$reference_dir
out_dir = config$output_dir
resource_dir = config$resource_dir
rdfiles = paste(rd_data_dir, dir(rd_data_dir), sep="/")



sample_id <- "EGAN00001367061"
chrom <- 1
# child_mrm <- read.table("/lustre/scratch113/projects/ddd/users/as33/MrM_page/results/15249082.ANXX.paired158.ba3ac4b206.cram.mrm.adm3.txt",header=T)
child_mrm <- read.table("/lustre/scratch113/projects/ddd/users/as33/MrM_page/results/14257319.ANXX.paired158.7bccc3500d.cram.mrm.adm3.txt",header=T)
dad_mrm <- read.table("/lustre/scratch113/projects/ddd/users/as33/MrM_page/results/15249070.ANXX.paired158.8bd2f156f8.cram.mrm.adm3.txt",header=T)
mum_mrm <- read.table("/lustre/scratch113/projects/ddd/users/as33/MrM_page/results/15249058.ANXX.paired158.fc17a04624.cram.mrm.adm3.txt",header=T)

child_chrom <- child_mrm[child_mrm$Chr == chrom,]
dad_chrom <- dad_mrm[dad_mrm$Chr == chrom,]
mum_chrom <- mum_mrm[mum_mrm$Chr == chrom,]
calls <- read.table("/nfs/users/nfs_a/as33/test_PAGE.txt",header=T)


