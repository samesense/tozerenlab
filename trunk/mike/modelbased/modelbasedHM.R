modelbasedHM<-function(out,groups,trueCI,nchains=25,max.iterations=NULL){

require(gplots)
require(SparseM)

#Plot a heatmap of the posterior pairwise distribution of expression data
#clustered with model-based clustering

#Inputs
#out[niterations,8] = values of scan number, log-likelihood, K, alpha, lambda,
#	eta, a, b from every iteration of model based clustering	
#groups[niterations,n] = posterior samples of the group membership function
#trueCI[n] = phenotype information
#nchains = number of parallel chains used in execution 
#max.iterations = used to discard a number of iterations

#Outputs
#writes heatmap plot to file modelbasedMap.pdf

#colors for colorbar
color<-c(rep('white',14),'blue','green','yellow','orange','red')
#generate a heatmap of the pairwise posterior distribution 

niterations<-dim(groups)[1]
nsamples<-dim(groups)[2]
niterClass<-niterations/nchains

#resort classes by nsamples/class
ns_class<-table(trueCI) #number of samples/class
classes<-ns_class@dimnames$trueCI #class labels in same order as ns_class
newCI_IDX<-rep(0,nsamples-15)
count<-0

for(i in 1:length(ns_class)){
	class_i<-match(max(ns_class),ns_class)
	classLabel<-classes[class_i]
	classIDX<-(1:nsamples)[trueCI==classLabel]
	
	newCI_IDX[(count+1):(count+max(ns_class))]<-classIDX
	count<-count+max(ns_class)
	ns_class<-ns_class[-class_i]
	classes<-classes[-class_i]
}	
trueCI<-trueCI[newCI_IDX]
groups<-groups[,newCI_IDX]

#generate pairwise posterior distribution
PD<-pairwisePD2(out,groups,nchains,max.iterations)$PD
PD<-as.matrix(PD)

classLabels<-unique(trueCI)
colbar<-rep(NA,length(trueCI))
sepvec<-rep(0,length(classLabels)-1)

#Heatmap of PD
for(i in 1:length(classLabels)){
	classIDX<-trueCI==classLabels[i]
	if(i<length(classLabels)){
		sepvec[i]<-max((1:nsamples)[classIDX])	
	}
	classPD<-PD[classIDX,classIDX]

	#order rows by dendrogram
	Rowv<-rowMeans(classPD,na.rm=T)
	hcr<-hclust(dist(classPD))
	ddr<-as.dendrogram(hcr)
	ddr<-reorder(ddr,Rowv)
	rowInd<-order.dendrogram(ddr)
	colInd<-rowInd

	#reorder rows
	temp<-PD[classIDX,]
	PD[classIDX,]<-temp[rowInd,]

	#reorder cols
	temp<-PD[,classIDX]
	PD[,classIDX]<-temp[,colInd]

	#colorbars for heatmap
	colbar[classIDX]<-color[i]
}
#heatmap
col=colorRampPalette(c('blue','yellow'),space='Lab')
col=col(10)

pdf(file='modelbasedMap.pdf')
ATHM(PD,TrueCI=unique(trueCI),Rowv=F,Colv=F,dendrogram='none',scale='none',trace='none',colsep=sepvec,rowsep=sepvec,ColSideColors=colbar,RowSideColors=colbar,key=T,col=col,labRow="",labCol="",density.info='none')
dev.off()
}	

ATHM<-function (x,TrueCI=NULL,SCIDX=NULL,SCIDX2=NULL,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, col = "heat.colors", colsep, 
    rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
    notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE), densadj = 0.25, main = NULL, xlab = NULL, 
    ylab = NULL, ...) 
{
    #############################################################	
    ### heatmap code adapted from heatmap.2 in gplots package ###	
    #############################################################
	
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) 
        if (missing(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    if (length(breaks) == 1) {
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
            length = breaks)
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    else if (is.character(col) && length(col) == 1) 
        col <- do.call(col, list(ncol))
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[] <- ifelse(x < min.breaks, min.breaks, x)
    x[] <- ifelse(x > max.breaks, max.breaks, x)
    lmat <- rbind(4:3, 2:1)
    lhei <- lwid <- c(keysize, 4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
           nr) 
          stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
	
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
	image(rbind(1:nr),col=RowSideColors[rowInd],axes=F)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
        cellnote <- t(cellnote)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = ncol(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        for (i in colInd) {
            if (!is.null(vline)) {
                vline.vals <- scale01(vline, min.scale, max.scale)
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        for (i in rowInd) {
            if (!is.null(hline)) {
                hline.vals <- scale01(hline, min.scale, max.scale)
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }else{
	 plot.new()
	x1<-rep(0.55,14)
	y1<-c(0.72,0.49,0.33,0.24,0.184,0.135,0.08)
	ciIDX1<-seq(2,14,2)
	
	text(x1,y1,cex=1.6,labels=TrueCI[ciIDX1])
	 }
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }else{ 
		plot.new()
	       x2<-c(0.08,0.41,0.60,0.72,0.79,0.835,0.90)
		 y2<-rep(0.22,length(x2))
             text(x2,y2,cex=1.6,labels=TrueCI[c(1,3,5,7,9,11,13)],srt=90)	
	}
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 2, 5, 1), cex = 0.75)
        if (symkey) {
            max.raw <- max(abs(x), na.rm = TRUE)
            min.raw <- -max.raw
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
	  z <- seq(min.raw, max.raw, length = length(col))
		
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("")
		#title("Color Key")
    }
    else plot.new()
	

    invisible(list(rowInd = rowInd, colInd = colInd,z=z,breaks=breaks,min.raw=min.raw,max.raw=max.raw))
}

pairwisePD2<-function(out,groups,nchains,max.iterations=NULL){

require(SparseM)

#Given a set of gibbs samplers run in parallel, find the most likely 
#partition from each chain and average them together

#Inputs
#out[niterations,8] = values of scan number, log-likelihood, K, alpha, lambda,
#	eta, a, b from every iteration of model based clustering	
#groups[niterations,n] = posterior samples of the group membership function
#trueCI[n] = phenotype information
#nchains = number of parallel chains used in execution 
#max.iterations = used to discard a number of iterations

#Outputs 
#List with three components 
#PD[n,n] = pairwise posterior distribution
#chainPD[nchains,n] = maximum likelihood partition from each run 
#chainLL = maximum likelihood from each run 

niterations<-dim(groups)[1]
nsamples<-dim(groups)[2]
niterChain<-niterations/nchains

if(is.null(max.iterations)){
	max.iterations<-niterChain
}else{
	max.iterations<-max.iterations
}

logL<-out[,2]
#find cluster indicator from each partition with the max likelihood 
chainIDX<-1:max.iterations
maxLL<-rep(0,length(nchains))
clustChain<-matrix(0,nrow=nchains,ncol=nsamples)

for(i in 1:nchains){
	tempLL<-logL[chainIDX]
	tempGroups<-groups[chainIDX,]
	maxLL[i]<-tempLL[tempLL %in% max(tempLL)]
	clustChain[i,]<-tempGroups[tempLL %in% max(tempLL),]
	chainIDX<-chainIDX+niterChain
}

#create pairwise posterior distribution from max likelihood 
#partitions from each chain
tempPwise<-as.matrix.csr(0,nsamples,nsamples)
for(i in 1:nrow(clustChain)){
	pwiseChain<-as.matrix.csr(outer(clustChain[i,],clustChain[i,],'=='))
	tempPwise<-tempPwise+pwiseChain
}
pwisePD<-tempPwise/nrow(clustChain)

#outliers
if(any(apply(x,2,function(x) all(x<=0.5)))) print(outliers) 
	
return(list(PD=pwisePD,chainPD=clustChain,chainLL=maxLL))
}


