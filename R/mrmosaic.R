wapplyp <- function(x, ws=window, by = NULL, FUN = NULL, ...) {
	FUN <- match.fun(FUN)
	if (is.null(by)) by <- ws
	lenX <- length(x)
	# get indices
	SEQ1 <- seq(1, lenX - ws + 1, by = by)
	# get points in windows
	SEQ2 <- lapply(SEQ1, function(x) x:(x + ws - 1))
	# Parallelised lapply 
	OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
	OUT <- base:::simplify2array(OUT, higher = TRUE)
	return(OUT)
}

distance_weighted_man_u <- function(indices,column="",cases=data.frame(),background=data.frame()) {
	subset <- cases[indices,]
	#Compute 
	middle_position<- subset$Position[length(subset$Position)/2+1]
	distances      <- abs(middle_position-subset$Position) + 1
	# Threshold for tricube distance
	h              <- 500000
	norm_distances <- distances/h
	# vector-based min comparisons
	# norm_distance is bounded by 1/500000 to 1
	norm_distances <- pmin(1,norm_distances)
	# tricube is bounded by 0 and 1
	tricube        <- (1-norm_distances^3)^3
	# get weight (proportion of all tricubed distances seen in this instance)
	w              <- tricube/sum(tricube)
	values         <- subset[,column]
	# higher prob of sampling values close to the center of the window
	cases_sample   <- sample(values,size=100,prob=w,replace=TRUE)
	bkg_sample     <- sample(background[,column],100,replace=TRUE)
	# H1: cases have a higher <column> than background
	man_wh         <- wilcox.test(cases_sample,bkg_sample,alternative="greater")$p.value
	return(man_wh)
}

get_mdpt_indices <- function (xdf, ws, by) {
	l    <- dim(xdf)[1]
	# start at mid pt of first window
	from <- 1+(ws-1)/2
	# to mid pt of last window
	to   <- (l - ws + 1 + (ws-1)/2 )
	is   <- seq(from=from, to=to, by=by)
	return (is)
} 

annotate_devs <- function(mydf=mydf,width=width,by=by){
        mydf <- na.omit(mydf)

        # logR smoothing for all genotypes
        mydf$logrdev        <- abs(mydf$Log.R.Ratio)
		# Here, rollapply is faster than multi-cored mclapply. #TODO: find alternative multi-threaded method
        mydf$rollmedianlogr <- rollapply(mydf$logrdev,width,median, fill=NA)
        
        # want the bdev smooth for the het sites
        mydf                <- subset(mydf,GType=="AB")
        mydf$bdev <- abs(0.5 - mydf[,"B.Allele.Freq"])
        return(mydf)
}

compute_stats  <- function(mydf=mydf,background=background,window=window,by=by,width=width,c=1) {
	
	# say("starting compute_stats, annotation")

	mydf       <- annotate_devs(mydf,      width=width, by=by)	
	background <- annotate_devs(background,width=width, by=by)

	indices    <- seq(1:dim(mydf)[1])
	is         <- get_mdpt_indices(mydf, ws=window, by=by)

	# say("finished annotation, starting first wp")

	manwh_rollmedianlogr <- wapplyp(indices,ws=window,FUN=distance_weighted_man_u,column="logrdev",cases=mydf,background=background,by=by)

	# say("finished first wp, starting second")

	manwh_bdev           <- wapplyp(indices,ws=window,FUN=distance_weighted_man_u,column="bdev",   cases=mydf,background=background,by=by)
	
	# say("finished second wp, starting manwh")

	mydf[is,"manwh_combined"] <- fisher.method(data.frame(manwh_bdev,manwh_rollmedianlogr))$S

	mydf <- na.omit(mydf)

	# say("finished manwh combined, returning form compute_stats")

	return(mydf)
}

merge_nearby_gada_segments <- function(mydf) {

	add_jcol_to_df <-function (mydf) {
		df2merge<-NULL; 
		j<-1; 
		df2merge<-rbind(data.frame(j=1,mydf[1,])); 
		for (i in 2:(dim(mydf)[1])) { 
			diff <- mydf[i,"IniPos"] - mydf[i-1,"EndPos"] ; 
			if (diff < 1e6) {
				df2merge <- rbind(df2merge,data.frame(j=j,mydf[i,])) 
			} else { 
				j<-j+1
				df2merge <- rbind(df2merge,data.frame(j=j,mydf[i,])) 
			} 
		}
		return(df2merge)
	}
	
	mymerge_base <- function(mydf,j) {
		a <- mydf[mydf$j==j,]
		chr   <- a[1,"Chrom"]
		start <- head(a,n=1)[,"IniPos"]
		end   <- tail(a,n=1)[,"EndPos"]
		lrr <- sum(a$LRR*a$nProbe)/sum(a$nProbe)
		# bdev only present in 'calling' not in 'simulations'
		if ("bdev" %in% colnames(a)) {
			bdev<- sum(a$bdev*a$nProbe)/sum(a$nProbe)
		}
		gada   <- a[1,"Tgada"]
		probe  <- sum(a$nProbe)
		gadaamp<- sum(a$GADAamp*a$nProbe)/sum(a$nProbe)
		if ("bdev" %in% colnames(a)) {
			data.frame(j=j,Chrom=chr,IniPos=start,EndPos=end,LRR=lrr,bdev=bdev,Tgada=gada,nProbe=probe,GADAamp=gadaamp)
		} else {
	
			data.frame(j=j,Chrom=chr,IniPos=start,EndPos=end,LRR=lrr,Tgada=gada,nProbe=probe,GADAamp=gadaamp)
		}
	}
	if (dim(mydf)[1] < 2) {
		a0 <- data.frame(mydf,Njoined=1)
		return(a0)
	}

	df2merge <- add_jcol_to_df(mydf)
	multi    <- table(df2merge$j)[which(table(df2merge$j)>1)]
	multi    <- data.frame(j=as.numeric(names(multi)),times=multi)
	a0 <- data.frame()
	if ( any(! df2merge$j %in% multi$j) ) a0 <- data.frame(df2merge[! df2merge$j %in% multi$j,],Njoined=1)
	if (dim(mydf)[1] > dim(a0)[1]) {
		for (i in 1:dim(multi)[1]) {
			j <- multi[i,"j"]
			Njoined <- multi[i,"times"]
			a0 <- rbind(a0 , data.frame(mymerge_base(df2merge,j),Njoined=Njoined) )
		}
	}
	a0 <- a0[order(a0$j),names(a0)!="j"]
	return(a0)

}

add_clon <- function(lrr=lrr,bdev=bdev) { 
    thresh <- 0.05
    if (lrr < -thresh)                  return(list(clon=loss_clon(bdev),type="loss"))
    if (lrr >  thresh)                  return(list(clon=gain_clon(bdev),type="gain"))
    if (lrr >= -thresh & lrr <= thresh) return(list(clon=loh_clon(bdev), type="loh"))
    stop()
}

loss_clon <- function(bdev) 2*bdev / (0.5 + bdev)
gain_clon <- function(bdev) 2*bdev / (0.5 - bdev)
loh_clon  <- function(bdev) 2*bdev

mrmosaic_calling_func <- function(data,t=t,window=40,width=5,by=1) {
    #Get the name of the sim file out of the path and use it for parsing out the coordinates
    CHROMOSOMES <- c(seq(1:22))
    output <- NULL	
    for (chr in CHROMOSOMES) {
    	print(paste("Processing chromosome ",chr))
        chrsubset <- subset(data,Chr==chr)
    	idx <- which(data$Chr != chr)
    	background <- data[sample(idx,10000),]
    	chrsubset <- compute_stats(mydf=chrsubset,background=background,window=window,by=by,width=width)
        #Run GADA segmentation on the Combined Fisher's Omnibus Statistic
        sink("/dev/null")
		Tgada  <- t
		dataSim<- setupGADAgeneral(chrsubset$manwh_combined)
		step1  <- SBL(dataSim, estim.sigma2=TRUE,maxit=1e7)
		step2  <- BackwardElimination(step1,T=Tgada,MinSegLen=15)
		s      <- summary(step2)
		events <- subset(s, State == 1)
		if(dim(events)[1] != 0){ #Check that there are 'significant' events being called otherwise skip
			events$IniPos <- chrsubset[events$IniProbe,]$Position
			events$EndPos <- chrsubset[events$EndProbe,]$Position
			events$Chrom  <- chr
			events$Tgada  <- t
			events$nProbe <- events$LenProbe
			events$GADAamp<- events$MeanAmp
			events$LRR    <- apply(cbind(events$Chrom, events$IniPos, events$EndPos), 
				MARGIN=1, FUN = function(x) { 
					mean(data$Log.R.Ratio[(data$Chr==x[1] & data$Pos >= x[2] & data$Pos <= x[3])]) 
			} )
			events$bdev   <- apply(cbind(events$IniProbe, events$EndProbe),
				MARGIN=1, FUN = function(x) {
					mean(  chrsubset$bdev[ x[1]:x[2] ]  )
			} )
			to_merge      <- events[,c("Chrom","IniPos","EndPos","LRR","bdev","Tgada","nProbe","GADAamp")]
			merged        <- merge_nearby_gada_segments(to_merge)
			tmp_return    <- apply(cbind(merged$LRR,merged$bdev), MARGIN=1, FUN = function(x) add_clon(lrr=x[1],bdev=x[2]) )
			tmp_output    <- cbind(merged,do.call(rbind, tmp_return))
			tmp_output$clon <- as.vector(unlist(tmp_output$clon))
			tmp_output$type <- as.vector(unlist(tmp_output$type))
			tmp_output2   <- tmp_output[,c(1,2,3,4,5,10,11,7,8)]
			output        <- rbind(output,tmp_output2)
			#output        <- rbind(output,merge_nearby_gada_segments(to_merge))
		}
       sink()
	}
	results <- list(data=data,sig.results=output)
	return(results)
}

annot_ll_and_recalc_sig <- function (results,resource_dir) {
	amp_file = paste(resource_dir,"ecdf_GADAamp.merged.gada.10.R3.Robj",sep="")
	probecount_file = paste(resource_dir,"ecdf_nprobes.merged.gada.10.R3.Robj",sep="")
	if(!file.exists(amp_file)){
		stop(paste("Could not find file (based on resource_dir in config): ", amp_file))
	}
	if(!file.exists(probecount_file)){
		stop(paste("Could not find file (based on resource_dir in config): ", probecount_file))
	}
	load(amp_file)
	load(probecount_file)
	
	tmp     <- results$sig.results
	
	tmp$p_nprobes <- 1 - ecdf_nprobes(tmp$nProbe)
	tmp$p_GADAamp <- 1 - ecdf_GADAamp(tmp$GADAamp)
	tmp$ll        <- log2(tmp$p_nprobes) + log2(tmp$p_GADAamp)
	
	results$sig.results <- tmp
	rm(tmp)
	return (results)
}
plot_lattice <- function (df_all, df_sig, sampleid, out_dir) {
	if (missing(df_all))   stop("Missing df_all")
	if (missing(df_sig))   stop("Missing df_sig")
	if (missing(sampleid)) stop("Missing sampleid")
	if (missing(out_dir))  stop("Missing out_dir")
#   if (is.na(list.dirs(out_dir))) system(cmd_to_make_dir)
    file_out     <-paste(sampleid,"all_chrs","png",sep=".")
    file_path_out<-paste(out_dir,file_out,sep="/")

    png(file=file_path_out,width=1200,height=700,res=75)

    par(mfrow=c( 3, 8 ), oma=(c(0,0,5,0)))
    chrs    <- c(seq(1,22),"X","Y")
    for (i in 1:length(chrs)) {
	    chr2plot     <-as.character(chrs[i])
		if ( length(which(df_all$chr==chr2plot)) != 0 ) {
			print(c("plotting lattice chromosome number",chr2plot))
			lattice_base(df_all, chr2plot)
			# df_sig header is:
			# Chrom   IniP    EndP    MeanAmp LR      Tgada 
	        df_sig_chr   <-df_sig[df_sig$Chrom==chr2plot,]
	        if (dim(df_sig_chr)[1] >= 1) plot_line_segs(df_sig_chr)
	    	title(paste("Chromosome", chr2plot),line=-1)
	    	title(outer=T,in_file_prefix,cex.main=5)
		} else { print(c("no chr data to plot for chromosome",chr2plot)) }
    }
    dev.off()
}

# Generic single-chromosome plotting function
# takes in a dataframe with positions to plot, 
# and a df of significant segments to plot on this chromosome
plot_single_chr_of_interest <- function (df_all, df_sig, sampleid, out_dir, chr2plot, ptsize=0.2, pch=19) {
	if (missing(df_all))   stop("Missing df_all")
	if (missing(df_sig))   stop("Missing df_sig")
	if (missing(sampleid)) stop("Missing sampleid")
	if (missing(out_dir))  stop("Missing out_dir")
	if (missing(chr2plot)) stop("Missing chr2plot")

    out_file     <-paste(sampleid,chr2plot,"png",sep=".")
    out_path_file<-paste(out_dir,out_file,sep="/")
	sampleid_chr <-paste(sampleid,chr2plot,sep=".")

	png(file=out_path_file,width=1200,height=700,res=75)
	single_chr_base_plot(df_all, chr2plot, pos, ptsize=ptsize,title=sampleid,pch=pch)
	df_sig_chr  <- df_sig[df_sig$Chrom==chr2plot,]
	if (dim(df_sig_chr)[1] >= 1) plot_line_segs(df_sig_chr)
    dev.off()
}
