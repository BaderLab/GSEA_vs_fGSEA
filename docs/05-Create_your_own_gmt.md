---
params:
  working_dir: ./generated_data/
  species: felis
  ensembl_dataset: fcatus_gene_ensembl
---

# Create GMT file from Ensembl

The [Baderlab geneset download site](https://download.baderlab.org/EM_Genesets/) is an updated resource for geneset files from GO, Reactome, WikiPathways, Pathbank, NetPath, HumanCyc, IOB, ... many others that can be used in [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) or [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) and many other enrichment tools that support the gmt format.  

Unfortunately genesets are only supplied for:

  * [Human](https://download.baderlab.org/EM_Genesets/current_release/Human/)
  * [Mouse](https://download.baderlab.org/EM_Genesets/current_release/Mouse/)
  * [Rat](https://download.baderlab.org/EM_Genesets/current_release/Rat/) 
  * [Woodchuck](https://download.baderlab.org/EM_Genesets/current_release/Woodchuck/)

If you are working in a different species you will need to generate your own gmt file. The best way to do this is through ensembl.  Ensembl doesn't have annotations for all the pathway databases listed above but it has annotations for most species from GO.


The parameters are set in the params option on this notebook but you can also manually set them here.

``` r
# for example - working_dir <- "./genereated_data"
working_dir <- params$working_dir

# for example - species <- "horse"
species <- params$species

# for example - ensembl_dataset <- "ecaballus_gene_ensembl"
ensembl_dataset <- params$ensembl_dataset
```



``` r
#use library
#make sure biocManager is installed
tryCatch(expr = { library("BiocManager")}, 
         error = function(e) { 
           install.packages("BiocManager")}, 
         finally = library("BiocManager"))


tryCatch(expr = { library("biomaRt")}, 
         error = function(e) { 
           BiocManager::install("biomaRt")}, 
         finally = library("biomaRt"))
```
## Load Libraries

Create or set a directory to store all the generatd results

``` r
if(!dir.exists(params$working_dir)){
  dir.create(params$working_dir)
}
```

## Set up Biomart connection

Connect to Biomart

``` r
ensembl <- useEnsembl(biomart = "genes", host = "asia.ensembl.org")
#ensembl <- useEnsembl("ensembl")
```


Figure out which dataset you want to use - for some species there might be a few datasets to choose from.  Not all of the datasets have common namesa associated with them.  For example, if you search for 'yeast' nothing will be returned but if you look for Saccharomyces or cerevisiae  you will be able to find it.


``` r
all_datasets <- listDatasets(ensembl)

#get all the datasets that match our species definition
all_datasets[grep(all_datasets$description,
                  pattern=species,
                  ignore.case = TRUE),]
```

```
##                dataset                 description         version
## 68 fcatus_gene_ensembl Cat genes (Felis_catus_9.0) Felis_catus_9.0
```

If you know the ensembl dataset that you want to use you can specify it in the parameters above or grab from the above table the dataset of the species that you are interested in. 


``` r
ensembl = useDataset(ensembl_dataset,mart=ensembl)
```

## Get species GO annotations

Get the GO annotations for our species

``` r
go_annotation <- getBM(attributes = c("external_gene_name",
                                      "ensembl_gene_id",
                                      "ensembl_transcript_id",
                                      "go_id", 
                                      "name_1006", 
                                      "namespace_1003",
                                      "go_linkage_type"), 
                       filters=list(biotype='protein_coding'), mart=ensembl);

#get just the go biological process subset
#####
# Get rid of this line if you want to include all of go and not just biological process
#####
go_annotation_bp <- go_annotation[which(
  go_annotation$namespace_1003 == "biological_process"),]

#compute the unique pathway sets
go_pathway_sets <- aggregate(go_annotation_bp[,1:5],
                             by = list(go_annotation_bp$go_id),
                             FUN = function(x){list(unique(x))})

#unlist the go descriptions
go_pathway_sets$name_1006 <- apply(go_pathway_sets,1,FUN=function(x){
   paste(gsub(unlist(x$name_1006),pattern= "\"",
              replacement = ""),collapse = "")})
```

There are two identifiers that you can choose from in the above table
 * external_symbols
 * ensembl_ids
 
 Each of these is stored as a list in the dataframe.  In order to convert it to the right format for the gmt file we need to convert the list to string of tab delimited strings.  (unfortunately there is no streaightforward way to write out a dataframe's column of lists.)

``` r
go_pathway_sets[1:3,"external_gene_name"]
```

```
## [[1]]
## [1] "MGME1"    "MPV17"    "AKT3"     "SLC25A36" "MEF2A"    "SLC25A33" "OPA1"    
## 
## [[2]]
## [1] ""      "MMP23"
## 
## [[3]]
##  [1] "TNP1"  "XNDC1" "SIRT1" "TDP1"  "APLF"  "ERCC6" "XRCC1" "APTX"  "LIG4" 
## [10] "ERCC8"
```

``` r
go_pathway_sets[1:3,"ensembl_gene_id"]
```

```
## [[1]]
## [1] "ENSFCAG00000023195" "ENSFCAG00000015150" "ENSFCAG00000003700"
## [4] "ENSFCAG00000035204" "ENSFCAG00000022495" "ENSFCAG00000022044"
## [7] "ENSFCAG00000000695"
## 
## [[2]]
## [1] "ENSFCAG00000034109" "ENSFCAG00000008341"
## 
## [[3]]
##  [1] "ENSFCAG00000036860" "ENSFCAG00000044053" "ENSFCAG00000012373"
##  [4] "ENSFCAG00000005982" "ENSFCAG00000022636" "ENSFCAG00000031322"
##  [7] "ENSFCAG00000005222" "ENSFCAG00000008728" "ENSFCAG00000031985"
## [10] "ENSFCAG00000026531"
```

## Format results into GMT file

Convert column of lists to a tab delimited string of gene names

``` r
go_pathway_sets$collapsed_genenames <- apply(go_pathway_sets,1,
                                             FUN=function(x){
   paste(gsub(unlist(x$external_gene_name),pattern= "\"",
              replacement = ""),collapse = "\t")
})
```


Convert column of lists to a tab delimited string of gene names

``` r
go_pathway_sets$collapsed_ensemblids <- apply(go_pathway_sets,1,
                                              FUN=function(x){
   paste(gsub(unlist(x$ensembl_gene_id),pattern= "\"",
              replacement = ""),collapse = "\t")
})
```

The format of the GMT file is described [https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29](here) and consists of rows with the following

  * Name
  * Description
  * tab delimited list of genes a part of this geneset
  
Write out the gmt file with genenames


``` r
gmt_file_genenames <- go_pathway_sets[,c("Group.1","name_1006",
                                         "collapsed_genenames")]
colnames(gmt_file_genenames)[1:2] <- c("name","description") 

gmt_file_genenames$name <- paste(gmt_file_genenames$description, gmt_file_genenames$name,sep="_")
gmt_file_genenames$name <- gsub(gmt_file_genenames$name,pattern = " ",replacement = "_")
gmt_file_genenames$name <- gsub(gmt_file_genenames$name,pattern = ":",replacement = "_")

gmt_genenames_filename <- file.path(params$working_dir, paste(species,ensembl_dataset,"GO_genesets_GN.gmt",sep = "_"))

write.table(x = gmt_file_genenames,file = gmt_genenames_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
```

Write out the gmt file with ensembl ids

``` r
gmt_file_ensemblids <- go_pathway_sets[,c("Group.1","name_1006",
                                          "collapsed_ensemblids")]
colnames(gmt_file_ensemblids)[1:2] <- c("name","description") 

gmt_file_ensemblids$name <- paste(gmt_file_ensemblids$description, gmt_file_ensemblids$name,sep="_")
gmt_file_ensemblids$name <- gsub(gmt_file_ensemblids$name,pattern = " ",replacement = "_")
gmt_file_ensemblids$name <- gsub(gmt_file_ensemblids$name,pattern = ":",replacement = "_")


gmt_ensemblids_filename <- file.path(params$working_dir, paste(species,ensembl_dataset,"GO_genesets_esemblids.gmt",sep = "_"))

write.table(x = gmt_file_ensemblids,file = gmt_ensemblids_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
```
