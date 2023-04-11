---
params:
  working_dir: ./generated_data/g_profiler
  data_dir: ./data
  genelist_file: Supplementary_Table1_Cancer_drivers.txt
  max_gs_size: 350
  min_gs_size: 3
  min_intersection: 3
  organism: hsapiens
---

# Run g:profiler from R

Detailed instructions on how to run [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) programmatically from R

The parameters are set in the params option on this notebook but you can also manually set them here.

```r
# for example - working_dir <- "./genereated_data"
working_dir <- params$working_dir

data_dir <- params$data_dir

# for example - species <- "horse"
genelist_file <- params$genelist_file

# max size of the genesets for example -  350
max_gs_size <- params$max_gs_size

# max size of the genesets for example - 3
min_gs_size <- params$min_gs_size

#min intersection between your genelist and the geneset - for example 3
min_intersection <- params$min_intersection

# organism parameter used for g:profiler.  First letter of first word in species name followed by 
# the second word
# for example - hsapiens
organism <- params$organism
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
## 1 query_1        TRUE 1.426353e-37      5653        121               103
## 2 query_1        TRUE 3.391992e-36      5882        121               103
## 3 query_1        TRUE 7.172333e-36      3097        121                81
## 4 query_1        TRUE 2.511953e-35      5724        121               101
## 5 query_1        TRUE 7.298888e-35      3540        121                84
##   precision     recall    term_id source
## 1 0.8512397 0.01822041 GO:0031323  GO:BP
## 2 0.8512397 0.01751105 GO:0080090  GO:BP
## 3 0.6694215 0.02615434 GO:0031325  GO:BP
## 4 0.8347107 0.01764500 GO:0051171  GO:BP
## 5 0.6942149 0.02372881 GO:0010604  GO:BP
##                                                term_name effective_domain_size
## 1               regulation of cellular metabolic process                 21128
## 2                regulation of primary metabolic process                 21128
## 3      positive regulation of cellular metabolic process                 21128
## 4      regulation of nitrogen compound metabolic process                 21128
## 5 positive regulation of macromolecule metabolic process                 21128
##   source_order                                        parents
## 1         7549             GO:0019222, GO:0044237, GO:0050794
## 2        18924                         GO:0019222, GO:0044238
## 3         7551 GO:0009893, GO:0031323, GO:0044237, GO:0048522
## 4        14399                         GO:0006807, GO:0019222
## 5         4369             GO:0009893, GO:0043170, GO:0060255
```

Filter the table to include just the columns that are required for the generic enrichment map file results [GEM](https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files). Restrict the results to just the ones that have at least min_gs_size and less than max_gs_size terms and  min_intersection size include only the term_id, term_name, p_value (and p_value again because the p_value is actually the corrected p-value.  The output file does not contain the nominal p_value.  For down stream analysis though it is expected to have both a p-value and a q-value so just duplicate the q-value as both p-value and q-value)


```r
# filer by params defined above
enrichment_results <- subset(enrichment_results,term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))
```

In order to create a proper Generic enrichment results file we will need a copy of the gmt file used by g:Profiler. (also to create an Enrichment map).

Download the gmt file used for this analysis from g:profiler


```r
#the link to the gmt file is static no matter what version
gprofiler_gmt_url <- "https://biit.cs.ut.ee/gprofiler//static/gprofiler_full_hsapiens.name.gmt"

#get version info gprofiler as the gmt file is always associated with a specific version of g:profiler
gprofiler_version <- get_version_info(organism=organism)

gprofiler_gmt_filename <- file.path(working_dir,
                                            paste("gprofiler_full", organism,
                                                  gprofiler_version$gprofiler_version,sep="_",
                                                  ".name.gmt"))

download.file(url = gprofiler_gmt_url, destfile = gprofiler_gmt_filename)

#load in the g:profiler geneset file
capt_output <- capture.output(genesets_gprofiler <- GSA.read.gmt(filename = gprofiler_gmt_filename))

names(genesets_gprofiler$genesets) <- genesets_gprofiler$geneset.names
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
  genes_collapsed <- unlist(lapply(genes,FUN=function(x){paste(x,collapse = ",")}))
  genes_collapsed_df <- data.frame(term_id = names(genes), genes = genes_collapsed,stringsAsFactors = FALSE)
  return(genes_collapsed_df)
}
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
if(nrow(enrichment_results) >0){    

           #add phenotype to the results
          enrichment_results <- cbind(enrichment_results,1)
          
          # Add the genes to the genesets
          subset_genesets <- genesets_gprofiler$genesets[which(genesets_gprofiler$geneset.names %in% enrichment_results$term_id)]
          
          genes <- getGenesetGenes(query_set, subset_genesets)
          
          enrichment_results <- merge(enrichment_results,genes,by.x=1, by.y=1)
          
          colnames(enrichment_results) <- c("name","description","p-value","q-value","phenotype","genes")
          
}

#output the enrichment map file
write.table(enrichment_results, file = file.path(working_dir, "gprofiler_GEM_using_gprof_gmt.txt"),
                                                 row.names = FALSE, col.names = TRUE,quote = FALSE)
```


## Run g:profiler with your own genesets

Download the latest [Bader lab genesets](https://download.baderlab.org/EM_Genesets/current_release/Human/)




## Upload the gmt file to gprofiler

In order to use your own genesets with g:Profiler you need to upload the the file to their server first.  The function will return an ID that you need to specify in the organism parameter of the g:Profiler gost function call. 

```r
custom_gmt <- upload_GMT_file(gmtfile=dest_gmt_file)
```

```
## Your custom annotations ID is gp__15KB_4Xzp_hzo
## You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
```

```
## Just use: gost(my_genes, organism = 'gp__15KB_4Xzp_hzo')
```

For this query we are specifying - 

  * query - the set of genes of interest, as loaded in from the Supplementary_Table1_Cancer_drivers.txt file.
  * significant - set to FALSE because we want g:Profiler to return all the results not just the ones that it deems significant by its perdetermined threshold.
  * ordered_query - set to TRUE because for this set of genes they are ordered in order of their significance
  * correction_method - set to fdr.  by default g:Profiler uses g:Scs
  * organism - set to the custom_gmt ID ( for this run it is - gp__15KB_4Xzp_hzo) that we received when we uploaded our genetset file.



```r
gprofiler_results_custom <- gost(query = query_set ,
                                     significant=FALSE,
                                 ordered_query = TRUE,
                                    exclude_iea=FALSE,
                                     correction_method = "fdr",
                                 organism = custom_gmt
                                     )
```

```
## Detected custom GMT source request
```



```r
 #get the gprofiler results table
enrichment_results_customgmt <- gprofiler_results_custom$result
    
enrichment_results_customgmt[1:5,]
```

```
##     query significant      p_value term_size query_size intersection_size
## 1 query_1        TRUE 6.730012e-37      4645        110                94
## 2 query_1        TRUE 2.579228e-35      4747        110                93
## 3 query_1        TRUE 2.579228e-35      4881        110                94
## 4 query_1        TRUE 4.280880e-34      5230        110                95
## 5 query_1        TRUE 2.069151e-33      3435        110                81
##   precision     recall
## 1 0.8545455 0.02023681
## 2 0.8454545 0.01959132
## 3 0.8545455 0.01925835
## 4 0.8636364 0.01816444
## 5 0.7363636 0.02358079
##                                                                          term_id
## 1                       REGULATION OF CELLULAR METABOLIC PROCESS%GOBP%GO:0031323
## 2              REGULATION OF NITROGEN COMPOUND METABOLIC PROCESS%GOBP%GO:0051171
## 3                        REGULATION OF PRIMARY METABOLIC PROCESS%GOBP%GO:0080090
## 4                  REGULATION OF MACROMOLECULE METABOLIC PROCESS%GOBP%GO:0060255
## 5 REGULATION OF NUCLEOBASE-CONTAINING COMPOUND METABOLIC PROCESS%GOBP%GO:0019219
##                                                  source
## 1 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol
## 2 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol
## 3 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol
## 4 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol
## 5 Human_GOBP_AllPathways_no_GO_iea_April_02_2023_symbol
##                                                        term_name
## 1                       regulation of cellular metabolic process
## 2              regulation of nitrogen compound metabolic process
## 3                        regulation of primary metabolic process
## 4                  regulation of macromolecule metabolic process
## 5 regulation of nucleobase-containing compound metabolic process
##   effective_domain_size source_order parents
## 1                 18525        18677    NULL
## 2                 18525        11634    NULL
## 3                 18525        11073    NULL
## 4                 18525        10657    NULL
## 5                 18525         5967    NULL
```

Filter the table to include just the columns that are required for the generic enrichment map file results [GEM](https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files). Restrict the results to just the ones that have at least min_gs_size and less than max_gs_size terms and  min_intersection size include only the term_id, term_name, p_value (and p_value again because the p_value is actually the corrected p-value.  The output file does not contain the nominal p_value.  For down stream analysis though it is expected to have both a p-value and a q-value so just duplicate the q-value as both p-value and q-value)


```r
# filer by params defined above
enrichment_results_customgmt <- subset(enrichment_results_customgmt,
                                       term_size >= min_gs_size & 
                                   term_size <= max_gs_size & 
                                   intersection_size >= min_intersection , 
                                 select = c(term_id,term_name,p_value,p_value ))
```


## Create enrichment Results files 

In order to use our results down stream in the Enrichment map we need to generate results files that we can pass to Enrichment Map.  

Load in the GMT file



## Create an output file of the results - Generic enrichment Map file from Baderlab gmt

The file requires - 

  * name
  * description
  * p-value
  * q-value
  * phenotyp
  * list of genes (overlap of query set and original geneset)
  
  The list of genes needs to be calculated using the gmt file and original query set.  For each geneset found in the result find the overlap between the set of genes that are a part of the geneset and the query set. 



```r
if(nrow(enrichment_results_customgmt) >0){    

           #add phenotype to the results
          enrichment_results_customgmt <- cbind(enrichment_results_customgmt,1)
          
          # Add the genes to the genesets
          subset_genesets <- genesets_baderlab_genesets$genesets[which(genesets_baderlab_genesets$geneset.names %in% enrichment_results_customgmt$term_id)]
          
          genes <- getGenesetGenes(query_set, subset_genesets)
          
          enrichment_results_customgmt <- merge(enrichment_results_customgmt,genes,by.x=1, by.y=1)
          
          colnames(enrichment_results_customgmt) <- c("name","description","p-value","q-value","phenotype","genes")
          
}

#output the enrichment map file
write.table(enrichment_results_customgmt, 
            file = file.path(working_dir, "gprofiler_GEM_using_gprof_gmt.txt"),
                                                 row.names = FALSE, col.names = TRUE,quote = FALSE)
```

