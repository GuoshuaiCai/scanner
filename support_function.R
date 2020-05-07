# ----------------------
# Supporting Functions
# ----------------------

library(future)
library(SIBER)


#subset data, the default cell number is 100
subset_data <- function(dataobj, ncell=500, method="ellipse") {
    clusters<-as.factor(dataobj@active.ident)
    cluster.vec<-levels(clusters)
    cell.sample<-NULL
    tsne<-as.matrix(Embeddings(dataobj, reduction = "tsne"))

    for(cluster in cluster.vec){
        cells<-colnames(dataobj)[which(clusters==cluster)]
        if(method=="sampling"){
            if(length(cells)>ncell){
                cells<-sample(cells,ncell)
            }
        }
        if(method=="ellipse"){
            p<-ncell/length(cells)
            if(p>1){
                p<-1
            }
            Y<-tsne[which(clusters==cluster),,drop=FALSE]
            mu <- colMeans(Y)
            Sigma <- cov(Y)
            Z <- pointsToEllipsoid(Y, Sigma, mu)
            inside <- ellipseInOut(Z, p = p)
            cells<-cells[inside]
        }
        cell.sample<-c(cell.sample,cells)
    }
    dataobj <- subset(dataobj,cells=cell.sample)
    dataobj
}

# calculate geneset enrichment
gs.activity <- function(dataobj, genelist, condition="All", method="Exp") {
    object_data   <- as.matrix(GetAssayData(dataobj))

    if(method=="Exp"){
        object_data   <- object_data
    }

    if(method=="Rnk"){
        object_data   <- apply(object_data,2,function(x) rank(x))
    }

    genes.checked <- intersect(genelist, rownames(dataobj))
    gs.cell <- NULL

    if(method=="EigGen"){      
        object_data <- calEigenGene(seuratobj=dataobj, genelist=genelist)
        genes.checked<-1
    }

    if(method=="ES"){
        object_data <- calES(seuratobj=dataobj, genelist=genelist)
        genes.checked<-1
    }

    if (length(genes.checked) > 0) {
        gs.cell <- t(object_data[genes.checked,, drop = FALSE])
        gs.cell <- as.data.frame(as.matrix(gs.cell))
    }

    gs.cell<-matrix(rowMeans(gs.cell),ncol=1,dimnames=list(rownames(gs.cell),"score"))

    if(condition!="All"){
        cluster<-paste(Idents(dataobj),dataobj@meta.data[[condition]],sep="_")
    }else{
        cluster<-Idents(dataobj)
    }

    gs.cluster<-list()
    for (i in unique(cluster)) {
        cells.i <- which(cluster==i)
        data.i <- gs.cell[cells.i,,drop=FALSE]
        if (length(cells.i) > 1) {
            data.i <- colMeans(data.i,na.rm=TRUE)
        }
        gs.cluster[[i]] <- data.i
    }
    gs.cluster<-do.call("rbind",gs.cluster)

    gs.exp<-list("cell"=gs.cell,"cluster"=gs.cluster)
    
    gs.exp
}


#EigenGene calculation
calEigenGene <- function(seuratobj=NULL, genelist=NULL) {
    obj.s<-as.matrix(GetAssayData(seuratobj))
    obj.s<-t(obj.s[which(rownames(obj.s) %in% genelist),,drop=FALSE])

    idx <- which(apply(obj.s,2,var)==0)
    if(length(idx)>0){
        obj.s<-obj.s[,-idx,drop=FALSE]
    }
    pc1<-prcomp(obj.s, center=TRUE, scale = TRUE)$x[,"PC1"]
    pc1<-matrix(pc1,nrow=1)
    colnames(pc1)<-names(seuratobj@active.ident)
    rownames(pc1)<-"EigenGene"
    return(pc1)
}

#ES calculation function
ES.fun<-function (genelist=NULL, gene.score=NULL, weighted.type = 0, type="expr",signed=FALSE){
    ngene <- length(gene.score)
    if(!signed){
        gene.score<-abs(gene.score)
    }
    if(type=="rank"){
        gene.score<-rank(gene.score)
    }
    sort.idx <- order(gene.score, decreasing = TRUE)
    names(sort.idx)<-names(gene.score)
    gene.score.sorted <- gene.score[sort.idx]
    cumsum.score <- matrix(0, nrow = ngene, ncol = 1)

    ngene.hit <- length(genelist)
    ngene.miss <- ngene - ngene.hit
    sort.idx.hit <- sort.idx[match(genelist, names(sort.idx))]
    cumsum.score[, 1] <- -1/ngene.miss
    if (weighted.type == 0) {
        cumsum.score[sort.idx.hit, 1] <- 1/ngene.hit
    }
    else if (weighted.type == 1) {
        gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit]
        cumsum.score[sort.idx.hit, 1] <- gene.score.sorted.hit/sum(gene.score.sorted.hit)
    }
    else {
        gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit]^weighted.type
        cumsum.score[sort.idx.hit, 1] <- gene.score.sorted.hit/sum(gene.score.sorted.hit)
    }
    cumsum.score <- apply(cumsum.score, 2, cumsum)
    t(cumsum.score)
}

#ES calcualtion
calES <- function(seuratobj=NULL, genelist=NULL) {
    gene.score<-as.matrix(GetAssayData(seuratobj))
    
    ES<-matrix(NA,nrow=1,ncol=ncol(gene.score))
    for(i in 1:ncol(gene.score)){
        score.i<-ES.fun(genelist=genelist, gene.score=gene.score[,i])
        ES[1,i]<-max(score.i)
    }
    #browser()
    #ES.scale <- t(x = scale(x = t(x = ES)))
    ES.scale <- t(x = t(x = ES))
    #ES.scale <- MinMax(data = ES.scale, max = 2.5, min = (-1) * 2.5)
    colnames(ES.scale)<-names(seuratobj@active.ident)
    rownames(ES.scale)<-"EnrichementScore"
    return(ES.scale)
}