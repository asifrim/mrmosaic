#!/usr/bin/env /software/R-3.3.0/bin/Rscript

VERSION=0.1

library(argparse)
library(yaml)
library(gada)
library(zoo)
library(mrmosaic)


####################################################################################################
## ARGUMENT PARSER


#Parser
parser <- argparse::ArgumentParser(description=paste("MrMosaic (version ",VERSION,"):"))
parser$add_argument("--input", help = "Path to input file (output of cnv_baf_data)")
parser$add_argument("--config",help= "Path to config file")


args <- parser$parse_args()

if (is.null(args$config) || is.null(args$input)){
	parser$print_help()
	stop("Error: You must provide an --input file and a --config file")
}

config <- yaml::yaml.load_file(args$config)
input_file = args$input
in_file = tail(strsplit(input_file, "/")[[1]],n=1)
rd_data_dir = config$reference_dir
out_dir = config$output_dir
resource_dir = config$resource_dir
rdfiles = paste(rd_data_dir, dir(rd_data_dir), sep="/")

####################################################################################################
## MAIN
print(paste(Sys.time(),": Calculating Log2 Ratios..."))
data = read_unique_depth_bait_data(input_file, length(count.fields(input_file)))
log2_ratio = calc_log2_ratio(data, rdfiles)
adm3_scores = calc_adm3_scores( data.frame( data, log2_ratio ) )
mrm_data = reformat_for_mrmosaic(input_file,adm3_scores)
output_path = paste(out_dir, paste(in_file,"adm3.txt",sep="."), sep="/")
write.table(mrm_data, file=output_path, sep="\t", row.names=F, quote=F)

print(paste(Sys.time(),": Running MrMosaic..."))
results <- mrmosaic_calling_func(mrm_data,t=10)
results <- annot_ll_and_recalc_sig(results,resource_dir)

in_file_prefix  = gsub(x=in_file,pattern=".extraction.txt",replacement="")
out_file_prefix = in_file_prefix
opf1 = paste(out_dir, paste(out_file_prefix, "results.txt",sep="."), sep="/")

write.table(results$sig.results, file=opf1, sep="\t", row.names=F, quote=F)

# Prepare second output file (ll_thresh)
df_all = data.frame(chr=results$data$Chr,pos=results$data$Pos,baf=results$data$B.All,lrr=results$data$Log.R)
df_sig = results$sig.results
df_sig = df_sig[df_sig$ll < -8 | df_sig$ll == "-Inf",]	# NOTE: Hardcoded
opf2   = paste(out_dir, paste(out_file_prefix, "sig_results.txt",sep="."), sep="/")
write.table(df_sig, file=opf2, sep="\t", row.names=F, quote=F)
print(paste(Sys.time(),": Finished..."))
