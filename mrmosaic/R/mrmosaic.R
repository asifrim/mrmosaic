
#' Applies a function using a rolling window approach
#'
#' @param x vector containing data
#' @param ws window size
#' @param by increment size (e.g. a value of 1 computes the function on a window around every point)
#' @param FUN the function to be applied on every window
#' @param ... arguments passed to FUN
#'
#' @return output vector
#' @export
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


#' Distance weighted Mann-Whitney U Test
#'
#' This function computes a Mann-Whitney U Test at each of the provided indices. The data is weighed based on
#' the genomic distance to the index to be tested (points closer to the position under investigation have higher weights).
#' Weighing is done based on the tricube distance (similar to loess local regression methods). Data points at distances higher than 0.5 Mb are considered considered to have a weight = 0.
#' The Mann-Whitney test is then performed using weighted resampling strategy.
#'
#' @param indices indices at which to compute the test
#' @param column index or name of column on which to compute the test
#' @param cases data.frame containing data for for case to be tested for deviations from background
#' @param background data.frame containing sample of the background distribution (e.g. data points on other chromosome/samples)
#'
#' @return Vector with one-sided p-values of distance w
#' @export
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


#' Get all mid-points of all windows for a vector of positions. Handles windows at the beginning and end of the vector.
#'
#' @param xdf vector of positions
#' @param ws window size
#' @param by increment size
#'
#' @return vector with mid-points of windows
#' @export
get_mdpt_indices <- function (xdf, ws, by) {
	l    <- dim(xdf)[1]
	# start at mid pt of first window
	from <- 1+(ws-1)/2
	# to mid pt of last window
	to   <- (l - ws + 1 + (ws-1)/2 )
	is   <- seq(from=from, to=to, by=by)
	return (is)
}

#' Computes the rolling median averaged Log2Ratio and B-Allele freq deviation. This method computes the rolling median average of the LogR deviation
#' in order to use information in homozygous sites to reduce variance.
#'
#' @param mydf Data.frame with Log2ratios and B-allele frequencies
#' @param width Width of the window to be used to compute the rolling median average of
#' @param by Increment size
#'
#' @return input dataframe with additional logrdev, rollmedianlogr and bdev columns
#' @export
annotate_devs <- function(mydf=mydf,width=width,by=by){
        mydf <- na.omit(mydf)

        # logR smoothing for all genotypes
        mydf$logrdev        <- abs(mydf$Log.R.Ratio)
		# Here, rollapply is faster than multi-cored mclapply. #TODO: find alternative multi-threaded method
        mydf$rollmedianlogr <- zoo::rollapply(mydf$logrdev,width,median, fill=NA)

        # want the bdev smooth for the het sites
        mydf                <- subset(mydf,GType=="AB")
        mydf$bdev <- abs(0.5 - mydf[,"B.Allele.Freq"])
        return(mydf)
}


#' Wrapper function to compute probability of significant copy-number and b-allele frequence deviations and combining those using Fisher's method.
#'
#' @param mydf data frame with data points of the sample to be tested
#' @param background data frame with data points representing the background distributions (e.g. other chromosome/control sample)
#' @param window window size for computing distance-weighted Mann-Whitney U tests for significance of deviations.
#' @param by increment size
#' @param width window size to compute rolling media average of logR deviations.
#'
#' @return input data.frame with additional "manwh_combined" column containing the test statistic (S) for the Fisher's combined method.
#' @export
compute_stats  <- function(mydf=mydf,background=background,window=window,by=by,width=width) {

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


#' Merges overlapping segments called by GADA
#'
#' @param mydf data.frame with GADA result segments.
#'
#' @return data.frame with merged segments
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


#' Returns the putative type of event and clonality for an event
#'
#' @param lrr logR ratio deviation
#' @param bdev  bdev deviation
#'
#' @return list with putative clonality and event type ("loh","gain","loss")
#' @export
#'
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


#' Super function that runs the whole Mr. Mosaic pipeline.
#'
#' @param data data.frame with the data (i.e. output of reformat_for_mrmosaic())
#' @param t T-value for the GADA segmentation algorithm
#' @param window window size for computing distance-weighted Mann-Whitney U tests for significance of deviations.
#' @param by increment size
#' @param width window size to compute rolling media average of logR deviations.
#'
#' @return list containing the used data and called structural variants
#' @export
mrmosaic_calling_func <- function(data,t=t,window=40,width=5,by=1) {
    #Get the name of the sim file out of the path and use it for parsing out the coordinates
    CHROMOSOMES <- c(seq(1:22))
    output <- NULL
    for (chr in CHROMOSOMES) {
    	print(paste("Processing chromosome ",chr))
        chrsubset <- subset(data,Chr==chr)
    	idx <- which(data$Chr != chr)
    	background <- data[sample(idx,10000,replace=TRUE),]
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


#' Annotate with loglikelihood of the event being a false positive given the number of probes encompassing the event and the amplitude of the GADA signal.
#'
#' @param results output of mr_mosaic_calling_func()
#' @param resource_dir path to directory containing the distributions of probecounts and amplitudes of false positive calls.
#'
#' @return input results file where sig.results has the probabilities of an event being a false positive given the number of probes and the amplitude. The joint loglikelihood is given in the "ll" column.
#' @export
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
