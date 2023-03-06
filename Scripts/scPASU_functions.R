## Adapted from https://github.com/hwanglab/apa_atingLab2019 ##

library(Rsamtools)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
library(BSgenome.Hsapiens.UCSC.hg38)
library(matrixStats)
library(edgeR)
library(DEXSeq)
library(GenomicRanges)
library(handy)
library(goldmine)
library(rtracklayer)
library(dplyr)


### This function creates TUs ###

reduceGenes <- function(genes,chrs,flank=5000,ncore=1)
{
  
	genes <- makeGRanges(genes,strand=T)
	genes <- genes[seqnames(genes) %in% chrs]
	genes$cdsStart <- genes$cdsStart+1
	genes$cdsEnd <- genes$cdsEnd+1
	stopifnot(genes$strand %in% c("+","-"))

	# Some genes may have disjoint txUnits, let's treat them as separate genes to make things easier
	message("Clustering TUs")
	df <- as.data.frame(table(genes$gene.id))
	genes.singles <- genes[genes$gene.id %in% df[df$Freq==1,]$Var1]
	genes.tored <- genes[genes$gene.id %in% df[df$Freq>1,]$Var1]
	spl <- split(genes.tored,genes.tored$gene.id)
	spl.red <- lapply(spl,function(x) reduce(x,with.revmap=FALSE))#collect union boundary of the same gene
	for(i in 1:length(spl.red))
	{
		spl.red[[i]]$gene.id <- names(spl)[i]
	}
	
	genes.red <- do.call(c,unname(spl.red))
	genes.red$name <- genes[match(genes.red$gene.id,genes$gene.id)]$name
	message("Single isoform annotated: ",length(unique(genes.singles$gene.id)))
	message("Multiple isoform annotated: ",length(unique(genes.tored$gene.id)))	
	dt <- makeDT(genes.red)
	genes.red.counts <- dt[,list(nPos=length(strand),nStrands=length(unique(strand))),by="gene.id"]
	message("Reduced to disjoint clusters on same strand only: ",nrow(genes.red.counts[(nPos>1)&(nStrands==1),]))
	message("Reduced to disjoint clusters on different strands: ",nrow(genes.red.counts[(nStrands>1),]))

	# Assign isoforms to TUs
	tu <- genes.singles
	values(tu) <- NULL # removes all metadata
	tu$gene.id <- genes.singles$gene.id
	tu$name <- genes.singles$name
	tu <- c(tu,genes.red)
	seqlevels(tu) <- chrs
	tu <- tu[order(tu)]
	tu$tu <- paste0("TU",1:length(tu))
    tu$tu_anno<-paste0(tu$tu,':',tu$name)

	# From now on, must aggregate by TU rather than gene.id, as there can be duplication in gene.id	
	message("Assigning transcripts to TUs")
	genes.singles$tu <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu
    genes.singles$tu_anno <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu_anno

	fo <- data.table(as.data.frame(findOverlaps(genes.tored,tu)))#to see how multiple-isoform genes are overlapped with the tu-assigned set; as a default, strand information is also included in the criteria
	fo$gene_gid <- genes.tored[fo$queryHits]$gene.id
	fo$tu_gid <- tu[fo$subjectHits]$gene.id
	fo$assign <- fo$gene_gid==fo$tu_gid
	fo <- fo[assign==TRUE,] #For the one w/ FALSE, a different gene is overlaped. The action with assign==TRUE still retains the one w/ overlapped with the other gene
	stopifnot(!duplicated(fo$queryHits))
	stopifnot(nrow(fo)==length(genes.tored))
	dt <- makeDT(genes.tored)
	dt[fo$queryHits,tu:= tu[fo$subjectHits]$tu] #bring tu id
    dt[,tu_anno:= paste0(dt$tu,':',dt$name)] #Add tu anno

	genes.tu <- c(genes.singles,makeGRanges(dt,strand=TRUE))

	stopifnot(length(unique(genes.tu$tu))==length(tu))
	stopifnot(length(genes.tu)==length(genes))
	stopifnot(length(unique(genes$gene.id))==length(unique(genes.tu$gene.id)))
	seqlevels(genes.tu) <- chrs
	genes.tu <- genes.tu[order(genes.tu)]

	# Assign TUs as coding
	message("Detecting coding TUs")
	genes.tu$coding <- genes.tu$cdsStart!=genes.tu$cdsEnd
	dt <- makeDT(genes.tu)
	dt <- dt[,list(coding=any(coding)),by="tu"]
	stopifnot(nrow(dt)==length(tu))
	tu$coding <- "new"
	tu[match(dt$tu,tu$tu)]$coding <- dt$coding

	# Get 3' end flanks (i.e., 5k bp segment right after TSE)
	message("Generating TU flanks")
	tu.flank <- flank(tu,width=flank,start=FALSE,both=FALSE,ignore.strand=FALSE)

	# Get 3' UTRs
	message("Clustering 3' UTRs for each coding TU")
	coding <- genes.tu[genes.tu$coding==TRUE]

	# If on (+), then want range between cdsEnd and txEnd
	# If on (-), then want range between txStart and cdsStart
	coding.p <- coding[strand(coding)=="+"]
	utr3.p <- GRanges(seqnames(coding.p),IRanges(coding.p$cdsEnd,end(coding.p)),strand=strand(coding.p),tu=coding.p$tu,gene.id=coding.p$gene.id,name=coding.p$name)
	coding.m <- coding[strand(coding)=="-"]
	utr3.m <- GRanges(seqnames(coding.m),IRanges(start(coding.m),coding.m$cdsStart),strand=strand(coding.m),tu=coding.m$tu,gene.id=coding.m$gene.id,name=coding.m$name)
	stopifnot((length(coding.p)+length(coding.m))==length(coding))

	redme <- function(myutrs)
	{			
		df <- as.data.frame(table(myutrs$tu))
		myutrs.singles <- myutrs[myutrs$tu %in% df[df$Freq==1,]$Var1]
		myutrs.tored <- myutrs[myutrs$tu %in% df[df$Freq>1,]$Var1]
		spl <- split(myutrs.tored,myutrs.tored$tu)
		spl.red <- lapply(spl,reduce)
		for(i in 1:length(spl.red))
		{
			spl.red[[i]]$tu <- names(spl)[i]
		}
		myutrs.red <- do.call(c,unname(spl.red))
		myutrs.red$gene.id <- myutrs[match(myutrs.red$tu,myutrs$tu)]$gene.id
		myutrs.red$name <- myutrs[match(myutrs.red$tu,myutrs$tu)]$name
		c(myutrs.singles,myutrs.red)
	}
	utr3.red.m <- redme(utr3.m)
	utr3.red.p <- redme(utr3.p)
	utr3 <- c(utr3.red.p,utr3.red.m)
	utr3 <- utr3[order(utr3)]
	utr3 <- utr3[width(utr3)>0]
	utr3$utr <- paste0("UTR",1:length(utr3))

	# Internal
	message("Computing coding gene internal regions via setdiff")	
	spl <- split(utr3,utr3$tu)
	coding <- tu[tu$coding==TRUE]
	
    # Remove coding genes that don't have 3' UTRs
	coding <- coding[coding$tu %in% names(spl)]
	diff <- lapply(coding$tu,function(x) setdiff(coding[coding$tu==x],spl[[x]]))
	for(i in 1:length(diff))
	{
		diff[[i]]$tu <- coding$tu[i]
	}
	
	internal <- do.call(c,diff)	
	internal$gene.id <- tu[match(internal$tu,tu$tu)]$gene.id
	internal$name <- tu[match(internal$tu,tu$tu)]$name

	# No-UTR gene bodies
	coding <- tu[tu$coding==TRUE]
	noutr <- coding[!(coding$tu %in% names(spl))]

	# Between annotated 3' UTRs
	message("Computing between multi-UTR")
	dt <- makeDT(utr3)
	dt <- dt[,list(chr=chr[1],start=min(start),end=max(end),strand=strand[1],gene.id=gene.id[1],name=name[1],nUtr=length(utr)),by="tu"]
	multis <- dt[nUtr>1,]
	# Now have the maximal range for all the multis
	# If we setdiff out the real 3', then we have just the between ranges
	union <- makeGRanges(multis,strand=T)
	spl <- split(utr3,utr3$tu)
	diff <- lapply(union$tu,function(x) setdiff(union[union$tu==x],spl[[x]]))
	for(i in 1:length(diff))
	{
		diff[[i]]$tu <- union$tu[i]
	}
	
	btw <- do.call(c,diff)	
	btw$gene.id <- tu[match(btw$tu,tu$tu)]$gene.id
	btw$name <- tu[match(btw$tu,tu$tu)]$name

	# Now remove the betweens from the internals
	internal <- internal[!(internal %in% btw)]
	
	ret <- list(tu=tu,genes.tu=genes.tu,flank=tu.flank,utr3=utr3,between=btw,internal=internal,noutr3=noutr)
	
    #writeBEDFromGRanges(ret$tu,file="output/tu_tu.bed",name="tu")
	#writeBEDFromGRanges(ret$flank,file="output/tu_flank.bed",name="tu")
	#writeBEDFromGRanges(ret$utr3,file="output/tu_utr3.bed",name="tu")
	#writeBEDFromGRanges(ret$between,file="output/tu_between.bed",name="tu")
	#writeBEDFromGRanges(ret$internal,file="output/tu_internal.bed",name="tu")
	#writeBEDFromGRanges(ret$noutr3,file="output/tu_noutr3.bed",name="tu")
	return(ret)
}

### This function assigns peaks to TUs ###

joinTus_peaks<- function(allpeaks,rg)
{
  # Retain only polya supported peaks
  keep<-which(allpeaks$polya=='polya')
  peaks<-allpeaks[keep,]

  #prs <- pr$pr

  # tu_red <- reduce(rg$tu)
  # fodt <- data.table(data.frame(findOverlaps(rg$tu,tu_red)))
  # fodt[,cnt:=.N,by=subjectHits]
  # rg$tu[fodt[cnt>1,queryHits],]
  # tu_red[fodt[cnt>1,subjectHits],]
  # tsv <- '/media/tommy/cache/exp_out/ating/polyASeqs/04_AnnoApa/output/overlapped_tu.tsv'
  # fwrite(file=tsv,data.table(data.frame(rg$tu[fodt[cnt>1,queryHits],])),sep="\t",quote=FALSE,row.names=FALSE)

        fo_tu <- data.table(as.data.frame(findOverlaps(peaks,rg$tu)))
        fo_flank <- data.table(as.data.frame(findOverlaps(peaks,rg$flank)))

        fo_tu$set <- "tu"
        fo_flank$set <- "flank"
        
        fo_tu$tu <- rg$tu[fo_tu$subjectHits]$tu
        fo_tu$coding <- rg$tu[fo_tu$subjectHits]$coding
        cat('Assigning TU anno \n')
	fo_tu$tu_anno <- rg$tu[fo_tu$subjectHits]$tu_anno
	cat('Printing overlap for tu \n')
	ft<-head(fo_tu)
	ft
	cat('\n')	


        fo_flank$tu <- rg$flank[fo_flank$subjectHits]$tu
        fo_flank$coding <- rg$flank[fo_flank$subjectHits]$coding
	cat('Assigning TU anno from flank \n')
        fo_flank$tu_anno <- rg$flank[fo_flank$subjectHits]$tu_anno
        cat('Printing overlap for flank \n')
        fo<-head(fo_flank)
        fo
	cat('\n')

     # Remove cases where same PR links to same TU and flank of that TU
        fo_flank$key <- paste0(fo_flank$queryHits,"+",fo_flank$tu)
        fo_tu$key <- paste0(fo_tu$queryHits,"+",fo_tu$tu)
        fo_flank <- fo_flank[!(fo_flank$key %in% fo_tu$key),]

        ov <- rbind(fo_tu,fo_flank)
	cat('Printing merged overlap: \n')
	o<-head(ov)
	o
	cat('\n')

        # Join table that links each PR to a TU/flank
        cat('Create data table \n')
	join <- data.table(peak=peaks[ov$queryHits]$peakID,tu=ov$tu,tu_anno=ov$tu_anno, type=ov$set,coding=ov$coding)

        # Count table for listing the contingencies
        agg <- join[,list(nTu=length(tu[type=="tu"]),nFlank=length(tu[type=="flank"])),by="peak"]

        # Reduced assignment table that assigns the uniques
        join <- join[,list(tu=tu,tu_anno=tu_anno,type=type,coding=coding,unique_peak=length(tu)==1,over_tus=toString(tu[type=="tu"]),flank_tus=toString(tu[type=="flank"])),by="peak"]
        
        join <- join[,list(peak=peak,type=type,coding=coding,unique_peak=unique_peak,unique_tu=all(unique_peak),over_tus=over_tus,flank_tus=flank_tus,tu_anno=tu_anno),by="tu"]

    summary <- data.frame(group="intergenic",n=length(peaks)-nrow(agg),stringsAsFactors=FALSE)
        summary <- rbind(summary,data.frame(group="unique_tu",n=sum((agg$nTu==1)&(agg$nFlank==0))))
        summary <- rbind(summary,data.frame(group="unique_flank",n=sum((agg$nTu==0)&(agg$nFlank==1))))
        summary <- rbind(summary,data.frame(group="multi_tu",n=sum((agg$nTu>1)&(agg$nFlank==0))))
        summary <- rbind(summary,data.frame(group="multi_flank",n=sum((agg$nTu==0)&(agg$nFlank>1))))
        summary <- rbind(summary,data.frame(group="unique_tu_multi_flank",n=sum((agg$nTu==1)&(agg$nFlank>0))))
        summary <- rbind(summary,data.frame(group="multi_tu_multi_flank",n=sum((agg$nTu>1)&(agg$nFlank>0))))
        summary$frac <- summary$n/sum(summary$n)
        stopifnot(sum(summary$n)==length(peaks))

        ret <-list(allpeaks=allpeaks,polya_peaks=peaks,join=join,agg=agg,summary=summary)
       return(ret)
}

#### Modified functions from () ###
		       
		       
getGenes2 <- function(geneset="ucsc",gencodetable=NULL,genome,cachedir=NULL,sync=TRUE,url=NULL)
{
    # Validate geneset
    if(!(geneset %in% c("ucsc","refseq","gencode","ensembl"))){stop("geneset must be one of \"ucsc\", \"refseq\", \"gencode\", or \"ensembl\"")}
    if((geneset=="gencode")&(is.null(gencodetable))){stop("Please provide the UCSC Genome Browser table name for the GENCODE gene build of interest (for example: w
gEncodeGencodeBasicV19)")}

    if(geneset=="ucsc")
    {
        kg <- getUCSCTable("knownGene",genome,cachedir,sync=sync)
        kgx <- suppressWarnings(getUCSCTable("kgXref",genome,cachedir,sync=sync))
        ki <- getUCSCTable("knownIsoforms",genome,cachedir,sync=sync)
        setnames(kg,"name","kgID")
        setnames(ki,"transcript","kgID")
        setkey(kg,kgID)
        setkey(kgx,kgID)
        setkey(ki,kgID)
        kg <- ki[kg,]
        genes <- kgx[kg,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=geneSymbol,gene.id=clusterId,isoform.id=NA,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds,kgID=kgID)]
        genes[,isoform.id:=kgID]
        genes[,kgID:=NULL]
        return(genes)
    } else if(geneset=="refseq")
    {
        rg <- getUCSCTable("refFlat",genome,cachedir,sync=sync)
        genes <- rg[,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=geneName,gene.id=geneName,isoform.id=name,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds)]
        return(genes)
    } else if(geneset=="ensembl")
    {
        eg <- getUCSCTable("ensGene",genome,cachedir,sync=sync)
        en <- getUCSCTable("ensemblToGeneName",genome,cachedir,sync=sync)
        setnames(eg,"name","isoform.id")
        setnames(en,"name","isoform.id")
        setkey(eg,"isoform.id")
        setkey(en,"isoform.id")
        genes <- en[eg,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=value,gene.id=name2,isoform.id=isoform.id,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds)]
      return(genes)
    } else if(geneset=="gencode")
    {
        gg <- getUCSCTable(gencodetable,genome,cachedir,sync=sync,url=url)
        genes <- gg[,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=name2,gene.id=name2,isoform.id=name,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds)]
    }
}



getTableHeaderFromSQL <- function(sql.file)
{
    cols <- c()

    con  <- file(sql.file, open = "r")

    extract <- FALSE
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
    {
        linevec <- unlist(str_split(oneLine," "))
        #print(linevec[3])
        #if(is.na(linevec[3])){linevec[3]<-"KEY"}
        #print(linevec[3])

        if((extract==TRUE)&((linevec[3]=="KEY")|(linevec[3]=="PRIMARY")|(linevec[3]=="UNIQUE")|(linevec[3]=="DEFAULT")))
        {
            extract <- FALSE
        }
        if(extract==TRUE)
        {
            #print(linevec[3])
            colname <- linevec[3]
            colname <- str_replace_all(linevec[3],"`","")
            cols <- c(cols, colname)
        }
        if(linevec[1]=="CREATE")
        {
            extract <- TRUE
        }
    }

    close(con)

    cols
}
