SCANNER takes Seurat object, which can be generated following the Seurat
pipleline. This instruction use a preprocessed Seurat object of data
from cells of human bronchial epithelium (Duclos et al. Sci Adv 2019).
Sequencing read counts in single cells were downloaded from NIH GEO
(GSE131391). Subsequent data analyses, including data normalization,
high variable feature selection, data scaling, dimension reduction, and
cluster identification, were performed using the Seurat 3.0 package.

Packages *Seurat*, *dplyr* and *SIBER* will be used

``` r
library("Seurat")
```

    ## Registered S3 method overwritten by 'R.oo':
    ##   method        from       
    ##   throw.default R.methodsS3

``` r
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library("SIBER")
```

Three files are required, including a Seurat object file
(GSE131391\_object.RData), an dataset description file
(GSE131391\_description.csv) and a phenotype file
(GSE131391\_series\_matrix.txt).

First, load the example Seurat object.

``` r
load("GSE131391_object.RData")
```

The dataset meta information should be prepared in the format same to
“GSE131391\_description.csv”. Include dataset meta information into
seurat object.

``` r
meta.info <- read.csv("GSE131391_description.csv",
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)
meta.info.l <- as.list(subset(meta.info, select = "Description", drop = T))
names(meta.info.l) <- row.names(meta.info)
seurat.object@misc$meta.info <- meta.info.l
```

Also, phenotype data can be obtained by reading the GEO series matrix
file or self-constructed matrix. Here, the information of age, sex and
smoking are obtained. Note that currently SCANNER only visualize groups
and thus two age groups (“&lt;30”,“&gt;=30”) were created.

``` r
sample.info <- GEOquery::getGEO(filename="GSE131391_series_matrix.txt")

age<-sample.info@phenoData@data$'age:ch1'
sex<-sample.info@phenoData@data$'Sex:ch1'
smoking<-sample.info@phenoData@data$'smoking status:ch1'
age.o<-which(as.numeric(age)>=30)
age.y<-which(as.numeric(age)<30)
age[age.o]<-">=30"
age[age.y]<-"<30"
```

Match and incoporte these phenotype information into the
*DataSegregation* of the Seurat object.

``` r
sample.vec<-sample.info@phenoData@data$title
match.loc<-match(gsub("_cell.+","",colnames(seurat.object)),gsub(" ","_",sample.vec))
age.vec<-factor(age[match.loc])
sex.vec<-factor(sex[match.loc])
smoking.vec<-factor(smoking[match.loc])
names(age.vec)<-colnames(seurat.object)
names(sex.vec)<-colnames(seurat.object)
names(smoking.vec)<-colnames(seurat.object)
seurat.object$age <-age.vec
seurat.object$sex <- sex.vec
seurat.object$smoking <- smoking.vec

seurat.object@misc$DataSegregation <- list(
            "age" =levels(seurat.object@meta.data$'age'),
            "sex" =levels(seurat.object@meta.data$'sex'),
            "smoking" =levels(seurat.object@meta.data$'smoking')
            )
```

Identify differentially expressed genes in each cell type compared to
others. Include the top 10, 30 and 100 into *misc* of the Seurat object.

``` r
DE.all <- FindAllMarkers(object = seurat.object,
                            only.pos = FALSE,
                            min.pct = 0.1)
```

    ## Calculating cluster club

    ## Calculating cluster goblet

    ## Calculating cluster basal

    ## Calculating cluster ciliated

    ## Calculating cluster WBC

    ## Calculating cluster basal-smoking

    ## Calculating cluster inocyte

``` r
seurat.object@misc$DE$top10<-DE.all %>% group_by(cluster) %>% 
                    top_n(10, avg_logFC) %>% as.data.frame()
seurat.object@misc$DE$top30<-DE.all %>% group_by(cluster) %>% 
                    top_n(30, avg_logFC) %>% as.data.frame()
seurat.object@misc$DE$top100<-DE.all %>% group_by(cluster) %>% 
                    top_n(100, avg_logFC) %>% as.data.frame()
```

**Optional but recommended** Subset data with a smaller number (500 by
default) of cells in each cell type. Two methods are aviable, including
(1) “sampling”: randomly select cells and (2) “ellipse”: select core
cells by fitting 2D ellipse. This subsetting step is automatally applied
with the “sampling” method when loading data into SCANNER. This
processing will lose cells but keep distinct ones. This step is
recommendated for largely reducing uploading and analysis time.

``` r
source("support_function.R")
```

    ## 
    ## Attaching package: 'future'

    ## The following object is masked from 'package:rmarkdown':
    ## 
    ##     run

``` r
seurat.object.subset<-subset_data(seurat.object,ncell=100,method="sampling")
```

Further, store the ranks of gene expression in each cell.

``` r
seurat.object.subset@assays$RNA@misc$rank<-
    apply(seurat.object.subset@assays$RNA@data,2,function(x) rank(x))
```

Last, save the data into a single R object (RDS format). It is ready to
load into SCANNER!

``` r
saveRDS(seurat.object.subset, "GSE131391_ellipse.RDS")
```

**Extra** We provide the functions for gene set activity analysis used
by SCANNER. They can be directly applied to the Seurat object. The gene
set activity in each cell can be inferred by four ways, including (1)
“Exp”: the average of gene expression, (2) “Rnk”: the average of Ranks
with lowest expression rank=1, (3) “EigGen”: the value of eigen gene and
(4) “ES”: the single sample gene enrichment score. The average activity
in each cluster (cell type) will be provided as well. Group comparison
can be easily achieved by letting *gs.activity* know which condition you
want to look at. In this example, we looked at the KEGG
Renin-angiotensin system pathway.

``` r
source("support_function.R")
ras<-c("ACE","ACE2","AGT","AGTR1","AGTR2","ANPEP","CMA1",
    "CPA3","CTSA","CTSG","ENPEP","LNPEP","MAS1","MME","NLN","REN","THOP1")
act.es<-gs.activity(seurat.object,ras,method="ES")
act.es$cluster
```

    ##                   score
    ## basal         0.5335443
    ## goblet        0.5296804
    ## ciliated      0.5304156
    ## basal-smoking 0.5324925
    ## WBC           0.5370387
    ## club          0.5299578
    ## inocyte       0.5300311

``` r
act.es.smoking<-gs.activity(seurat.object,ras,condition="smoking",method="ES")
act.es.smoking$cluster
```

    ##                                  score
    ## basal_Current Smoker         0.5371362
    ## goblet_Current Smoker        0.5294119
    ## ciliated_Current Smoker      0.5284983
    ## basal-smoking_Current Smoker 0.5328937
    ## WBC_Current Smoker           0.5362448
    ## club_Current Smoker          0.5291890
    ## inocyte_Current Smoker       0.5344121
    ## basal_Never Smoker           0.5318248
    ## ciliated_Never Smoker        0.5321151
    ## club_Never Smoker            0.5299792
    ## goblet_Never Smoker          0.5381845
    ## WBC_Never Smoker             0.5376130
    ## basal-smoking_Never Smoker   0.5188524
    ## inocyte_Never Smoker         0.5256502
