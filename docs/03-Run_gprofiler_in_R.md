---
params:
  working_dir: ./generated_data/g_profiler
  data_dir: ./data
  genelist_file: Supplementary_Table1_Cancer_drivers.txt
  max_gs_size: 250
  min_gs_size: 3
  min_intersection: 3
  organism: hsapiens
---

# Run g:profiler from R

## Initialize variables and libraries

Detailed instructions on how to run [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) programmatically from R

The parameters are set manually here but if you want to run the script from the command line then you can update the notebook to pull the parameters from the command line given arguments by updating each variable below to pull the values from the paramters - for example:

  * variable <- params$parameter_name
  
For more details see - [defining and using parameters](https://bookdown.org/yihui/rmarkdown/params-declare.html) and [Knitting with parameters](https://bookdown.org/yihui/rmarkdown/params-knit.html)



```r
#where to put all the generated files
working_dir <- "./generated_data/g_profiler"

# where to find the data files needed to run the analysis
data_dir <-  "./data"

# File name containing the list of genes to be used for analysis
genelist_file <- "Supplementary_Table1_Cancer_drivers.txt"

# default max size of the genesets for example -  250.  For this example we
# will be varying this parameter
max_gs_size <- 250

# default min size of the genesets for example -  3
min_gs_size <- 3

#min intersection between your genelist and the geneset - for example 3
min_intersection <- 3

# organism parameter used for g:profiler.  
# First letter of first word in species name followed by 
# the second word for example - hsapiens
organism <- "hsapiens"
```



```r
#use library
tryCatch(expr = { library("gprofiler2")}, 
         error = function(e) { 
           install.packages("gprofiler2")}, 
         finally = library("gprofiler2"))

tryCatch(expr = { library("GSA")}, 
         error = function(e) { 
           install.packages("GSA")}, 
         finally = library("GSA"))
```

Create or set a directory to store all the generatd results

```r
if(!dir.exists(params$working_dir)){
  dir.create(params$working_dir)
}
```

## Load in Query set

Load in the set of genes that we will be running g:profiler with

```r
 #load in the file
    current_genelist <- read.table(file = 
                                     file.path(data_dir, genelist_file),
                                   header = FALSE,
                                   sep = "\t", quote = "",
                                   stringsAsFactors = FALSE)

  query_set <- current_genelist$V1
```


With regards to pathway sets there are two options when using [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) - 

  * Use the genesets that are supplied by [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost)
  * Upload your own genesets. 
  
The most common reasons for supplying your own genesets is the ability to use up to date annotations or in-house annotations that might not be available in the public sphere yet.  One of the greatest features of [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) is that it is updated on a regular basis and most of the previous versions are available online ont the [gprofiler archive](https://biit.cs.ut.ee/gprofiler/page/archives).

The [gprofielr2](https://biit.cs.ut.ee/gprofiler/page/r) -[g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) R implementation is a wrapper for the web version.  You require an internet connection to get enrichment results.  

## Run g:profiler with supplied genesets

For detailed descriptions of all the parameters that can be specified for the gost g:profiler function see -[here](https://rdrr.io/cran/gprofiler2/man/gost.html)

For this query we are specifying - 

  * query - the set of genes of interest, as loaded in from the Supplementary_Table1_Cancer_drivers.txt file.
  * significant - set to FALSE because we want g:Profiler to return all the results not just the ones that it deems significant by its perdetermined threshold.
  * ordered_query - set to TRUE because for this set of genes they are ordered in order of their significance
  * correction_method - set to fdr.  by default g:Profiler uses g:Scs
  * organism - set to "hsapiens" for homo sapiens.  Organism names are constructed by concatenating the first letter of the name and the family name (according to gprofiler2 documentation)
  * source - the geneset source databases to use for the analysis.  We recommend using GO biological process (GO:BP), WikiPathways (WP) and Reactome (Reac) but there are additional sources you can add (GO molecular function or cellular component(GO:MF, GO:CC), KEGG, transcription factors (TF), microRNA targets (MIRNA), corum complexes (CORUM), Human protein atlas (HPA),Human phenotype ontology (HP) ) 


```r
gprofiler_results <- gost(query = query_set ,
                          significant=FALSE,
                          ordered_query = TRUE,
                          exclude_iea=FALSE,
                          correction_method = "fdr",
                          organism = organism,
                          source = c("REAC","WP","GO:BP"))
```



```r
 #get the gprofiler results table
enrichment_results <- gprofiler_results$result
    
enrichment_results[1:5,]
```

```
##     query significant      p_value term_size query_size intersection_size
## 1 query_1        TRUE 9.272751e-39      5493        121               103
## 2 query_1        TRUE 1.694623e-36      5836        121               103
## 3 query_1        TRUE 1.276992e-35      5662        121               101
## 4 query_1        TRUE 4.128404e-35      2976        121                79
## 5 query_1        TRUE 1.198524e-34      3131        121                80
##   precision     recall    term_id source
## 1 0.8512397 0.01875114 GO:0031323  GO:BP
## 2 0.8512397 0.01764907 GO:0080090  GO:BP
## 3 0.8347107 0.01783822 GO:0051171  GO:BP
## 4 0.6528926 0.02654570 GO:0031325  GO:BP
## 5 0.6611570 0.02555094 GO:0051173  GO:BP
##                                                    term_name
## 1                   regulation of cellular metabolic process
## 2                    regulation of primary metabolic process
## 3          regulation of nitrogen compound metabolic process
## 4          positive regulation of cellular metabolic process
## 5 positive regulation of nitrogen compound metabolic process
##   effective_domain_size source_order
## 1                 21110         7510
## 2                 21110        18790
## 3                 21110        14316
## 4                 21110         7512
## 5                 21110        14318
##                                          parents
## 1             GO:0019222, GO:0044237, GO:0050794
## 2                         GO:0019222, GO:0044238
## 3                         GO:0006807, GO:0019222
## 4 GO:0009893, GO:0031323, GO:0044237, GO:0048522
## 5             GO:0006807, GO:0009893, GO:0051171
```

## Download and load g:profiler geneset file

In order to create a proper Generic enrichment results file we will need a copy of the gmt file used by g:Profiler. (also to create an Enrichment map).

Download the gmt file used for this analysis from g:profiler


```r
#the link to the gmt file is static no matter what version
gprofiler_gmt_url <- 
  "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_hsapiens.name.gmt"

#get version info gprofiler as the gmt file is always associated with 
# a specific version of g:profiler
gprofiler_version <- get_version_info(organism=organism)

gprofiler_gmt_filename <- file.path(working_dir,
                                  paste("gprofiler_full", organism,
                                    gprofiler_version$gprofiler_version,sep="_",
                                    ".name.gmt"))

if(!file.exists(gprofiler_gmt_filename)){
  download.file(url = gprofiler_gmt_url, 
              destfile = gprofiler_gmt_filename)
}
```

To create a proper Generic enrichmentMap results file we need to include the list of genes that are associated with each geneset.  To do that we need to know what genes are associated with each set and filter them by our query set.  Load in the geneset definitions from the gmt file we just downloaded from g:profiler site.  


```r
#load in the g:profiler geneset file
capt_output <- capture.output(genesets_gprofiler <- GSA.read.gmt(
                                      filename = gprofiler_gmt_filename))

names(genesets_gprofiler$genesets) <- genesets_gprofiler$geneset.names
```

For the next module the name of the gmt file is - gprofiler_full_hsapiens.name.gmt but it is important to preserve the database version so in the future when we revisit these results for publication or results verfication we have the exact version used.  Instead of creating a copy of the file (which can be pretty large) create a symbolic link to the file with the generic name.


```r
#file.exists does not work for a symbolic link on my computer for some reason
# list the files in the directory and check if the symbolic link is there
#if(file.exists(file.path(working_dir, "gprofiler_full_hsapiens.name.gmt"))){
if(length(grep(x = list.files(file.path(working_dir)), 
              pattern = "gprofiler_full_hsapiens.name.gmt",
              fixed = TRUE) > 0 )){

  file.remove(file.path(working_dir, "gprofiler_full_hsapiens.name.gmt"))
}
```

```
## [1] TRUE
```

```r
file.symlink( gprofiler_gmt_filename,file.path(working_dir, 
                                   "gprofiler_full_hsapiens.name.gmt"))
```

```
## [1] TRUE
```



```r
# Given:
# query_genes - genes used for enrichment analysis (or as query)
#
# returns - the genes that overlap with the query set and part of the given
#           genesets
getGenesetGenes <- function(query_genes, subset_genesets){
  genes <- lapply(subset_genesets,FUN=function(x){intersect(x,query_genes)})
  
  # For each of the genes collapse to the comma separate text
  genes_collapsed <- unlist(lapply(genes,FUN=function(x){
                                                paste(x,collapse = ",")}))
  
  genes_collapsed_df <- data.frame(term_id = names(genes), 
                            genes = genes_collapsed,stringsAsFactors = FALSE)
  
  return(genes_collapsed_df)
}
```

## Filter results by geneset size

Filter the table to include just the columns that are required for the generic enrichment map file results [GEM](https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files). Restrict the results to just the ones that have at least min_gs_size and less than max_gs_size terms and  min_intersection size include only the term_id, term_name, p_value (and p_value again because the p_value is actually the corrected p-value.  The output file does not contain the nominal p_value.  For down stream analysis though it is expected to have both a p-value and a q-value so just duplicate the q-value as both p-value and q-value)

Vary the thresholds for max_gs_size just as we did in Module 2 lab - 

  * min_gs_size = 3
  
  * max_gs_size = 10000
  * max_gs_size = 1000
  * max_gs_size = 250
  


```r
# filer by params defined above
# by default we have set the max and min gs size to 250 and 3, respectively.
enrichment_results_mxgssize_250_min_3 <- 
                        subset(enrichment_results,term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))

enrichment_results_mxgssize_1000_min_3 <- 
                        subset(enrichment_results,term_size >= min_gs_size & 
                                   term_size <= 1000 & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))

enrichment_results_mxgssize_10000_min_3 <- 
                        subset(enrichment_results,term_size >= min_gs_size & 
                                   term_size <= 10000 & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))
```


## Create an output file of the results - Generic enrichment Map file from g:profiler gmt

The file requires - 

  * name
  * description
  * p-value
  * q-value
  * phenotyp
  * list of genes (overlap of query set and original geneset)
  
  The list of genes needs to be calculated using the gmt file and original query set.  For each geneset found in the result find the overlap between the set of genes that are a part of the geneset and the query set. 



```r
# Given:
# gprofiler_results - results form g_profiler R function (filtered by desired)
# parameters
# gs - genes associated with each geneset, loaded in from a gmt file. 
#
# returns - the properly formatted GEM file results
#
createGEMformat <- function(results, gs, query_genes){

  if(nrow(results) >0){    

           #add phenotype to the results
          formatted_results <- cbind(results,1)
          
          # Add the genes to the genesets
          subset_genesets <- gs$genesets[
            which(gs$geneset.names 
                  %in% results$term_id)]
          
          genes <- getGenesetGenes(query_genes, subset_genesets)
          
          formatted_results <- merge(formatted_results,genes,by.x=1, by.y=1)
          
          colnames(formatted_results) <- c("name","description","p-value",
                                           "q-value","phenotype","genes")
          
  }
  return(formatted_results)
}
```



```r
enrichment_results_mxgssize_10000_min_3_GEMfile <- createGEMformat(
  enrichment_results_mxgssize_10000_min_3, genesets_gprofiler, query_set)

enrichment_results_mxgssize_1000_min_3_GEMfile <- createGEMformat(
  enrichment_results_mxgssize_1000_min_3, genesets_gprofiler, query_set)

enrichment_results_mxgssize_250_min_3_GEMfile <- createGEMformat(
  enrichment_results_mxgssize_250_min_3, genesets_gprofiler, query_set)
```



Output each of the above filtered files


```r
#output the enrichment map file
write.table(enrichment_results_mxgssize_10000_min_3_GEMfile, 
            file = file.path(working_dir, 
                "gProfiler_hsapiens_lab2_results_GEM_maxterm10000.txt"),
            row.names = FALSE, 
            col.names = TRUE,
            quote = FALSE)

#output the enrichment map file
write.table(enrichment_results_mxgssize_1000_min_3_GEMfile, 
            file = file.path(working_dir, 
                "gProfiler_hsapiens_lab2_results_GEM_maxterm1000.txt"),
            row.names = FALSE, 
            col.names = TRUE,
            quote = FALSE)

#output the enrichment map file
write.table(enrichment_results_mxgssize_250_min_3_GEMfile, 
            file = file.path(working_dir, 
                "gProfiler_hsapiens_lab2_results_GEM_maxterm250.txt"),
            row.names = FALSE, 
            col.names = TRUE,
            quote = FALSE)
```




## Run g:profiler with your own genesets (example using BaderLab genesets)

## Download and load Bader lab geneset file

Download the latest [Bader lab genesets](https://download.baderlab.org/EM_Genesets/current_release/Human/)


```r
gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"

#list all the files on the server
filenames = RCurl::getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include 
# terms inferred from electronic annotations(IEA)
#start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
  contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(working_dir,gmt_file)

if(!file.exists(dest_gmt_file)){
  download.file(
    paste(gmt_url,gmt_file,sep=""),
    destfile=dest_gmt_file
  )
}
```

In order to use our results down stream in the Enrichment map we need to generate results files that we can pass to Enrichment Map.  

Load in the GMT file


## Filter Bader lab geneset file

The g:Profiler interface only allows for filtering genesets by size only after the analysis is complete.  After the analysis is complete means the filtering is happening after Multiple hypothesis testing.  Filtering prior to the analysis will generate more robust results because we exclude the uninformative large genesets prior to testing changing the sets that multiple hypothesis filtering will get rid of.  

Create multiple gmt files with different filtering thresholds - remove 
  * genesets greater than 250 genes
  * geneset greater than 1000 genes
  * geneset greater than 10000 genes
  

```r
# Filter geneset GSA object by specified gs size threshold 
#
# Given - 
# genesets - in GSA object
# gs_sizes - list of all the sizes of the genesets found in the genesets
# filter_threshold - value to filter the geneset by.  
# 
# returns - filtered genesets in GSA object
filter_genesets <- function(genesets, gs_sizes, filter_threshold) {
  
  filtered_genesets <- genesets
  
  filtered_genesets$genesets <- genesets$genesets[
                      which(gs_sizes<filter_threshold)]
  filtered_genesets$geneset.names <- genesets$geneset.names[
                      which(gs_sizes<filter_threshold)]
  filtered_genesets$geneset.descriptions <- genesets$geneset.descriptions[
                      which(gs_sizes<filter_threshold)]

  return(filtered_genesets)
}

# You can not simply write a list of lists to a file in R.  In order
# to output the new geneset file you need to convert it ot a data.frame
# To do this convert the list of genes to a tab delmiated list in one column
# of the dataframe.
# format to write out to a file. 
#
# Given - 
# genesets - in GSA object
 
# returns - formatted genesets as data frame
  
format_genesets <- function(genesets) {
    
  collapsed_genesets <- data.frame(name=genesets$geneset.names, 
                            description= genesets$geneset.description)
  collapsed_genesets$genes <- unlist(lapply(genesets$genesets,
                                             FUN=function(x){
                                              paste(x,collapse = "\t")
                                            }))
  
  return(collapsed_genesets)
}
```



The format of the GMT file is described [https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29](here) and consists of rows with the following

  * Name
  * Description
  * tab delimited list of genes a part of this geneset
  
Write out the gmt file with genenames


```r
#get the geneset sizes
gs_sizes_baderlab_sets <- lapply(genesets_baderlab_genesets$genesets,
                                 FUN = function(x){
                                               length(x)
                                                  })

# max 10,000
genesets_baderlab_genesets_max10000 <- filter_genesets(genesets_baderlab_genesets,
                                                gs_sizes_baderlab_sets, 
                                                10000)

genesets_baderlab_genesets_max10000_filename <- gsub(x =dest_gmt_file, 
                                                  pattern = "symbol" ,
                                                  replacement = "symbol_max10000"
                                                     )

if(!file.exists(genesets_baderlab_genesets_max10000_filename)){

  write.table(x = format_genesets(genesets_baderlab_genesets_max10000),
            file = genesets_baderlab_genesets_max10000_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
}

#max gs size of 1,000
genesets_baderlab_genesets_max1000 <- filter_genesets(genesets_baderlab_genesets,
                                                      gs_sizes_baderlab_sets, 
                                                      1000)
genesets_baderlab_genesets_max1000_filename <- gsub(x =dest_gmt_file, 
                                                pattern = "symbol" ,
                                                replacement = "symbol_max1000"
                                                     )
if(!file.exists(genesets_baderlab_genesets_max1000_filename)){

  write.table(x = format_genesets(genesets_baderlab_genesets_max1000),
            file = genesets_baderlab_genesets_max1000_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
}

#max gs size of 250
genesets_baderlab_genesets_max250 <- filter_genesets(genesets_baderlab_genesets,
                                                      gs_sizes_baderlab_sets, 
                                                      250)


genesets_baderlab_genesets_max250_filename <- gsub(x =dest_gmt_file, 
                                                  pattern = "symbol" ,
                                                  replacement = "symbol_max250"
                                                     )
if(!file.exists(genesets_baderlab_genesets_max250_filename)){
  write.table(x = format_genesets(genesets_baderlab_genesets_max250),
            file = genesets_baderlab_genesets_max250_filename,
            quote = FALSE,sep = "\t",row.names = FALSE,
            col.names=TRUE)
}
```


## Upload the gmt files to gprofiler

In order to use your own genesets with g:Profiler you need to upload the the file to their server first.  The function will return an ID that you need to specify in the organism parameter of the g:Profiler gost function call. 

```r
custom_gmt_max250 <- upload_GMT_file(
                        gmtfile=genesets_baderlab_genesets_max250_filename)
```

```
## Your custom annotations ID is gp__RBhF_ElIY_R28
## You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
```

```
## Just use: gost(my_genes, organism = 'gp__RBhF_ElIY_R28')
```

```r
custom_gmt_max1000 <- upload_GMT_file(
                        gmtfile=genesets_baderlab_genesets_max1000_filename)
```

```
## Your custom annotations ID is gp__Te0O_xWNH_aK8
## You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
```

```
## Just use: gost(my_genes, organism = 'gp__Te0O_xWNH_aK8')
```

```r
custom_gmt_max10000 <- upload_GMT_file(
                        gmtfile=genesets_baderlab_genesets_max10000_filename)
```

```
## Your custom annotations ID is gp__4ZB8_XKvw_Zug
## You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
```

```
## Just use: gost(my_genes, organism = 'gp__4ZB8_XKvw_Zug')
```

For this query we are specifying - 

  * query - the set of genes of interest, as loaded in from the Supplementary_Table1_Cancer_drivers.txt file.
  * significant - set to FALSE because we want g:Profiler to return all the results not just the ones that it deems significant by its perdetermined threshold.
  * ordered_query - set to TRUE because for this set of genes they are ordered in order of their significance
  * correction_method - set to fdr.  by default g:Profiler uses g:Scs
  * organism - set to the custom_gmt ID ( for this run it is - gp__RBhF_ElIY_R28) that we received when we uploaded our genetset file.



```r
gprofiler_results_custom_max250 <- gost(query = query_set ,
                                     significant=FALSE,
                                 ordered_query = TRUE,
                                    exclude_iea=FALSE,
                                     correction_method = "fdr",
                                 organism = custom_gmt_max250
                                     )
```

```
## Detected custom GMT source request
```

```r
gprofiler_results_custom_max1000 <- gost(query = query_set ,
                                     significant=FALSE,
                                 ordered_query = TRUE,
                                    exclude_iea=FALSE,
                                     correction_method = "fdr",
                                 organism = custom_gmt_max1000
                                     )
```

```
## Detected custom GMT source request
```

```r
gprofiler_results_custom_max10000 <- gost(query = query_set ,
                                     significant=FALSE,
                                 ordered_query = TRUE,
                                    exclude_iea=FALSE,
                                     correction_method = "fdr",
                                 organism = custom_gmt_max10000
                                     )
```

```
## Detected custom GMT source request
```



```r
 #get the gprofiler results table
enrichment_results_customgmt_max250 <- gprofiler_results_custom_max250$result
enrichment_results_customgmt_max1000 <- gprofiler_results_custom_max1000$result
enrichment_results_customgmt_max10000 <- gprofiler_results_custom_max10000$result

    
enrichment_results_customgmt_max250[1:5,]
```

```
##     query significant      p_value term_size query_size intersection_size
## 1 query_1        TRUE 1.385009e-22        64         78                17
## 2 query_1        TRUE 1.988040e-20        68         78                16
## 3 query_1        TRUE 6.846277e-15        54        101                13
## 4 query_1        TRUE 3.849668e-14        97        107                15
## 5 query_1        TRUE 1.546752e-13       159        108                17
##   precision    recall
## 1 0.2179487 0.2656250
## 2 0.2051282 0.2352941
## 3 0.1287129 0.2407407
## 4 0.1401869 0.1546392
## 5 0.1574074 0.1069182
##                                                                                   term_id
## 1         HEAD AND NECK SQUAMOUS CELL CARCINOMA%WIKIPATHWAYS_20220510%WP4674%HOMO SAPIENS
## 2               GLIOBLASTOMA SIGNALING PATHWAYS%WIKIPATHWAYS_20220510%WP2261%HOMO SAPIENS
## 3 PATHWAYS AFFECTED IN ADENOID CYSTIC CARCINOMA%WIKIPATHWAYS_20220510%WP3651%HOMO SAPIENS
## 4                                     CELL CYCLE%WIKIPATHWAYS_20220510%WP179%HOMO SAPIENS
## 5                          REGULATION OF CELL CYCLE G1/S PHASE TRANSITION%GOBP%GO:1902806
##                                                         source
## 1 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol_max250
## 2 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol_max250
## 3 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol_max250
## 4 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol_max250
## 5 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol_max250
##                                        term_name effective_domain_size
## 1          Head and neck squamous cell carcinoma                 17047
## 2                Glioblastoma signaling pathways                 17047
## 3  Pathways affected in adenoid cystic carcinoma                 17047
## 4                                     Cell cycle                 17047
## 5 regulation of cell cycle G1/S phase transition                 17047
##   source_order parents
## 1         5002    NULL
## 2         5554    NULL
## 3         4968    NULL
## 4         4934    NULL
## 5        19009    NULL
```

Filter the table to include just the columns that are required for the generic enrichment map file results [GEM](https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files). Restrict the results to just the ones that have at least min_gs_size and less than max_gs_size terms and  min_intersection size include only the term_id, term_name, p_value (and p_value again because the p_value is actually the corrected p-value.  The output file does not contain the nominal p_value.  For down stream analysis though it is expected to have both a p-value and a q-value so just duplicate the q-value as both p-value and q-value)


```r
# filer by params defined above
enrichment_results_customgmt_max250 <- subset(enrichment_results_customgmt_max250,
                                       term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))

enrichment_results_customgmt_max1000 <- subset(enrichment_results_customgmt_max1000,
                                       term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))

enrichment_results_customgmt_max10000 <- subset(enrichment_results_customgmt_max10000,
                                       term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))
```


## Create an output file of the results - Generic enrichment Map file from Baderlab gmt

Use the same function defined above but instead of passing the genesets from the g_profiler gmt file pass the geneset defitnions we loaded in from the Baderlab gmt file. 


```r
enrichment_results_customgmt_GEM_max250 <- createGEMformat(
                                    enrichment_results_customgmt_max250, 
                                    genesets_baderlab_genesets_max250, 
                                    query_set)

#output the enrichment map file
write.table(enrichment_results_customgmt_GEM_max250, 
                  file = file.path(
                    working_dir, "gProfiler_hspaiens_Baderlab_max250.txt"),
                  row.names = FALSE, 
                  col.names = TRUE,
                  quote = FALSE)
       
enrichment_results_customgmt_GEM_max1000 <- createGEMformat(
                                  enrichment_results_customgmt_max1000,
                                  genesets_baderlab_genesets_max1000, 
                                  query_set)

#output the enrichment map file
write.table(enrichment_results_customgmt_GEM_max1000, 
                  file = file.path(
                    working_dir, "gProfiler_hspaiens_Baderlab_max1000.txt"),
                  row.names = FALSE, 
                  col.names = TRUE,
                  quote = FALSE)

enrichment_results_customgmt_GEM_max10000 <- createGEMformat(
                                  enrichment_results_customgmt_max10000,
                                  genesets_baderlab_genesets_max10000, 
                                  query_set)

#output the enrichment map file
write.table(enrichment_results_customgmt_GEM_max10000, 
                  file = file.path(
                    working_dir, "gProfiler_hspaiens_Baderlab_max1000.txt"),
                  row.names = FALSE, 
                  col.names = TRUE,
                  quote = FALSE)
```

