####################################################################################################
# Generate log2 ratio and adm3 values from exome sequence read count data
# Author: Tomas William Fitzgerald
# Email: tf2@sanger.ac.uk

# reads output from cnv_baf_data.c uniqifed on cnv bait regions 

read_unique_depth_bait_data <- function(file, n) {
command = paste("awk -F", '"\t"', " '!_[$1]++' ", file, " | cut -f1-3 ",sep="")
read.table(pipe(paste("awk -F", '"\t"', " '!_[$1]++' ", file, " | cut -f1-3 ",sep="")), header=FALSE, check.names=FALSE, colClasses=c("character", "numeric", "numeric"), nrows=n, comment.char="", sep="\t")
}

# calculates the correlation between samples using read counts at bait regions
generate_correlation_values <- function(data, rdfiles) {
sapply(rdfiles, function(x) cor(data[,2]+0.01, read_unique_depth_bait_data(x, nrow(data))[,2]+0.01))
}

# generate log2 ratio values for a sample
calc_log2_ratio <- function(data, rdfiles) {
	cors = generate_correlation_values(data, rdfiles)
	auto_files = rdfiles[order(cors, decreasing=T)]
	auto_files=auto_files[2:(length(rdfiles))]
	big_data = data.frame(NA,do.call('cbind',lapply(auto_files, function(x) read_unique_depth_bait_data(x, n=nrow(data) ) [,2])))[,-(1)]+0.01
	med_data = apply(big_data, 1, median)
return(log2((data[,2]+0.01)/med_data))
}

# converts read counts and weights to adm3 scores
adm3_score <- function(d) {
	w = 1/(d[,5]^2);
	EWS = d[,4]*w/sqrt(w);
return(data.frame(d,w,EWS));
}

# splits the read count distribution into bins
breakpoints <- function(x,max_bin_size=1000) {
	xxc = data.frame(x[,1],as.numeric(x[,2])+0.01,as.numeric(x[,3])+0.01)
	max <- max_bin_size
	d1 <- split(sort(xxc[,2]), ceiling(seq_along(xxc[,2])/max))
	nmax <- c(1:length(d1))
	for(j in 1:length(d1)) { nmax[j] <- max(d1[[j]])  }
	nnmax <- c(0,nmax[1:length(nmax)-1],max(xxc[,2])+1)
	nnmax <- unique(nnmax)
	xc <- cut(xxc[,2],breaks=nnmax)
	xxc_cut <- data.frame(xxc,xc)
	xcu <- sort(unique(xc))	
	bpoints <- data.frame(xcu,nnmax[1:length(nnmax)-1]); names(bpoints) = c("#ReadsInRanges","#Breakpoints");
return(bpoints)
}

# generates weights based on mad of lg2r within read count range and calculate adm3 scores
calc_adm3_scores <- function(readdata) {
	bps = breakpoints(readdata)
	nnmax = c(bps[,2],max(readdata[,2])+1)
	xc = cut(readdata[,2]+0.01,breaks=nnmax)
	xxc_cut = data.frame(seq(1:length(readdata[,1])),readdata,xc)
	names(xxc_cut)[length(xxc_cut)] = "cr"
	xcu = sort(unique(xc))
	CR_MADs = data.frame(xxc_cut[,1:5],xxc_cut$cr)
	names(CR_MADs) = c("rowcount","position","rc","rd","lg2r","cr")	
	xclmad = rep(0.0001,length(xcu))
	for(i in 1:length(xcu)) {
		tmpc = CR_MADs[CR_MADs[,6]==xcu[i],]	
		xclmad[i] = mad(tmpc[,5],constant=1, na.rm=T)	
	}	
	BPranges = data.frame(xclmad,xcu)
	names(BPranges) = c("mads","cr")
	xtm = merge(x=CR_MADs,y=BPranges,by="cr")
	xtm = xtm[with(xtm,order(xtm[,2])),]
	adm_rep = adm3_score(xtm[,c(3:7)]);
	colnames(adm_rep) = c("position", "rc", "rd", "lg2r", "mads", "weight", "adm3")
return(adm_rep);
}