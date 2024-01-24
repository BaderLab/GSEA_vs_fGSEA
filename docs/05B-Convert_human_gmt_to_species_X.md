---
title: "Translate gmt file from human to another species"
output: html_notebook
params:
  working_dir: /home/rstudio/projects/generated_data/
  dest_gmt_file: 
  new_species: felis
  ensembl_dataset: fcatus_gene_ensembl
---

# Translate human gmt file to another species using ensembl

for species where there is a not an existing effort to characterize genes through GO, a much richer gmt file will come from mapping the gmt file from a well defined species preferrably the species your experiment is trying to model (for example, if you are modelling a human disease in mouse using a human gmt might return much richer results but species specific changes might be lost )

## Load in required libraries


```r
#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  
           install.packages("RCurl")}, 
         finally = library("RCurl"))

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

## Download the latest pathway definition file

Only Human, Mouse, Rat, and Woodchuck gene set files are currently available on the baderlab downloads site.  If you are working with a species other than human (and it is either rat,mouse or woodchuck) change the gmt_url below to the correct species. Check [here](http://download.baderlab.org/EM_Genesets/current_release/) to see all available species.

To create your own GMT file using Ensembl see [Create GMT file from Ensembl]


```r
dest_gmt_file <- "" #params$dest_gmt_file 

if(dest_gmt_file == ""){
  gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"
  
  #list all the files on the server
  filenames = getURL(gmt_url)
  tc = textConnection(filenames)
  contents = readLines(tc)
  close(tc)
  
  #get the gmt that has all the pathways and does not include terms 
  # inferred from electronic annotations(IEA)
  #start with gmt file that has pathways only and GO Biological Process only.
  rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
    contents, perl = TRUE)
  gmt_file = unlist(regmatches(contents, rx))
  
  dest_gmt_file <- file.path(params$working_dir,gmt_file )
  
  #check if this gmt file already exists
  if(!file.exists(dest_gmt_file)){
    download.file(
      paste(gmt_url,gmt_file,sep=""),
      destfile=dest_gmt_file
    )
  }
}
```



Load in the originaing gmt file


## Set up Biomart connection

Connect to Biomart

```r
ensembl <- useEnsembl(biomart = "genes" , host = "asia.ensembl.org")
#ensembl <- useEnsembl("ensembl")
```


Figure out which dataset you want to use - for some species there might be a few datasets to choose from.  Not all of the datasets have common namesa associated with them.  For example, if you search for 'yeast' nothing will be returned but if you look for Saccharomyces or cerevisiae  you will be able to find it.


```r
all_datasets <- listDatasets(ensembl, verbose=TRUE)
```

```
## Attempting web service request:
## asia.ensembl.org:80/biomart/martservice?redirect=no&type=datasets&requestid=biomaRt&mart=ENSEMBL_MART_ENSEMBL
```

```r
#get all the datasets that match our species definition
all_datasets[grep(all_datasets$description,
                  pattern=paste(params$new_species,sep=""),
                  ignore.case = TRUE),]
```

```
##                dataset                 description         version
## 68 fcatus_gene_ensembl Cat genes (Felis_catus_9.0) Felis_catus_9.0
```

Based on the above table define the dataset

```r
#get all the datasets that match our species definition
base_dataset <- all_datasets$dataset[grep(all_datasets$description,
                  pattern=paste(params$new_species,sep=""),
                  ignore.case = TRUE)]
```



If you know the ensembl dataset that you want to use you can specify it in the parameters above or grab from the above table the dataset of the species that you are interested in. 


```r
ensembl = useDataset(base_dataset,mart=ensembl)
```

## Convert the species genes from human

Get the homologs in the species of interest of the human genes.

```r
homologs <- getBM(attributes = c("external_gene_name",
                                      "ensembl_gene_id",
                                      "hsapiens_homolog_ensembl_gene", 
                                      "hsapiens_homolog_associated_gene_name" 
                                      ), 
                       filters=list(biotype='protein_coding'), mart=ensembl);
# get rid of rows that have an empty gene name for new species or for human
homologs <- homologs[which((homologs$external_gene_name != "") & (homologs$hsapiens_homolog_associated_gene_name != "")),]
```

Convert the human gmt file into the species of interest


```r
names(originating_gmt$genesets) <- originating_gmt$geneset.names

temp_genesets <- lapply(originating_gmt$genesets,FUN=function(x){homologs$external_gene_name[which(homologs$hsapiens_homolog_associated_gene_name %in% unlist(x))]})

translated_genessets <- list(geneset_names=originating_gmt$geneset.names, 
                             genesets = temp_genesets,
                             geneset.descriptions=originating_gmt$geneset.descriptions)
```


## Format results into GMT file

Convert lists to a tab delimited string of gene names

```r
translated_genessets$collapsed_genenames  <- unlist(lapply(translated_genessets$genesets,
                                             FUN=function(x){
   paste(unlist(x),collapse = "\t")
}))
```


The format of the GMT file is described [https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29](here) and consists of rows with the following

  * Name
  * Description
  * tab delimited list of genes a part of this geneset
  
  Write out the gmt file with genenames


```r
gmt_file_genenames <- data.frame(name=translated_genessets$geneset_names,
                                 description=translated_genessets$geneset.descriptions,
                                 genes = translated_genessets$collapsed_genenames)

gmt_genenames_filename <- file.path(params$working_dir, paste(params$new_species,"convertedfrom",basename(dest_gmt_file),sep = "_"))

write.table(x = gmt_file_genenames,file = gmt_genenames_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
```
