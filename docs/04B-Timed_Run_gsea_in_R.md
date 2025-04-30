---
params:
  analysis_name: BRCA_hd_tep
  working_dir: ./data/
  output_dir: ./generated_data/gsea/
  rnk_file: brca_hd_tep_ranks.rnk
  gsea_jar: /home/rstudio/GSEA_4.3.3/gsea-cli.sh
  gsea_directory: 'data/Human_GOBP_AllPathways_noPFOCR_no_GO_iea_May_01_2024_symbol.gmt'
  run_gsea: true
---
# Run GSEA from within R

This notebook is based largely on the [original notebook](https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/Protocol2_createEM.html) published with EnrichmentMap Protocol[@em2019] 

There is no package to run the original algorithm of GSEA[@gsea2005] in R.  There are many packages that have been published to imitate the process but none are recognized by The GSEA team.  
![GSEA R package message](./images/gsea_r_package_message.png)


## Load in required libraries


``` r
#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  
           install.packages("RCurl")}, 
         finally = library("RCurl"))


tryCatch(expr = { library("benchmarkme")}, 
         error = function(e) {  
           install.packages("benchmarkme")}, 
         finally = library("benchmarkme"))
```



## Configurable Parameters

In order to run GSEA automatically through the notebook you will need to download the gsea jar from [here](http://software.broadinstitute.org/gsea/downloads.jsp).  Specify the exact path to the gsea jar in the parameters in order to automatically compute enrichments using GSEA.

If you are running this notebook using the [baderlab workshop docker image](https://hub.docker.com/r/risserlin/workshop_base_image) then the image comes pre-installed with the gsea jar that you can use to run gsea directly in the docker.  The path to the GSEA jar in the docker is - /home/rstudio/GSEA_4.3.2/gsea-cli.sh

In order to run GSEA automatically you need to speciry the path to the gsea jar file.
The gsea_jar needs to be the full path to the GSEA 4.3.3 directory that you downloaded from GSEA. for example  /Users/johnsmith/GSEA_4.3.3/gsea-cli.sh

The parameters are set manually here but if you want to run the script from the command line then you can update the notebook to pull the parameters from the command line given arguments by updating each variable below to pull the values from the paramters - for example:

  * variable <- params$parameter_name
  
For more details see - [defining and using parameters](https://bookdown.org/yihui/rmarkdown/params-declare.html) and [Knitting with parameters](https://bookdown.org/yihui/rmarkdown/params-knit.html)


``` r
#path to GSEA jar 
# defined in the paramters at top of notebook
gsea_jar <- params$gsea_jar
```

Set the working directory as the directory to the directory where you downloaded all protocol files.  For example /User/JohnSmith/EMProtocolFiles/data


``` r
# defined in the paramters at top of notebook

#directory where all the data files are found.  For example -   ./data/ 
working_dir <- params$working_dir

#directory where all the data files are found.  For example -   ./generated_data/gsea/
output_dir <- params$output_dir
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

#The name to give the analysis in GSEA - for example Basal_vs_Classical
analysis_name <- params$analysis_name

#rank file to use in GSEA analysis.  
#For example - TCGA-PAAD_GDC_Subtype_Moffitt_BasalvsClassical_ranks.rnk
rnk_file <- params$rnk_file

#run_gsea - true/false
# This parameter is for the compilation of the notebook.  
run_gsea <- params$run_gsea

#set the gmt file you want to use if you don't want to use the latest gmt file.
# For example, if you set dest_gmt_file =="" the below script will automatically
# download the latest gmt file from baderlab webstie.  If it is set then it
# will use the file specified.  
dest_gmt_file = params$gsea_directory
```




## Download the latest pathway definition file

Only Human, Mouse, Rat, and Woodchuck gene set files are currently available on the baderlab downloads site.  If you are working with a species other than human (and it is either rat,mouse or woodchuck) change the gmt_url below to the correct species. Check [here](http://download.baderlab.org/EM_Genesets/current_release/) to see all available species.

To create your own GMT file using Ensembl see [Create GMT file from Ensembl]


``` r
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
  rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_noPFOCR_no_GO_iea.*.)(.gmt)(?=\">)",
    contents, perl = TRUE)
  gmt_file = unlist(regmatches(contents, rx))
  
  dest_gmt_file <- file.path(output_dir,gmt_file )
  
  #check if this gmt file already exists
  if(!file.exists(dest_gmt_file)){
    download.file(
      paste(gmt_url,gmt_file,sep=""),
      destfile=dest_gmt_file
    )
  }
} else {
  file.copy(dest_gmt_file,to = output_dir)
}
```

```
## [1] TRUE
```

***
## Run GSEA
(GSEA)[http://software.broadinstitute.org/gsea/index.jsp] is a stand alone java program with many customizable options.  It can be easily run through its integrated user interface.  To make this a seemless pipeline we can run GSEA from the command line with a set of options.  Any of the supplied options can be customized and there are many additional options that can be specified.  For more details see (here)[http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_GSEA_from]

In the below command the following options have been specified:

 * rnk - path to the rank file
 * gmx - path to the gene set definition (gmt) file
 * collapse - true/false indicates whether the expression/rnk file needs to be collapsed from probes to gene symbols
 * nperm - number of permutations
 * scoring_scheme - 
 * rpt_label - name of the directory with output
 * rnd_seed - random seed to use
 * set_max - maximum size for individual gene sets.  In GSEA interface this is set to 500 but we prefer to use a more stringent setting of 200. 
 * set_min - minimum size for individual gene sets 
 * zip_report - true/false to zip output directory
 * out - directory where to place the result directory.

 

``` r
start_time <- Sys.time()

if(run_gsea){
  command <- paste("",gsea_jar,  
                   "GSEAPreRanked -gmx", dest_gmt_file, 
                   "-rnk" ,file.path(working_dir,rnk_file), 
                   "-collapse false -nperm 1000 -scoring_scheme weighted", 
                   "-rpt_label ",analysis_name,
                   "  -plot_top_x 20 -rnd_seed 12345  -set_max 500",  
                   " -set_min 15 -zip_report false ",
                   " -out" ,output_dir, 
                   " > gsea_output.txt",sep=" ")
  system(command)
}

end_time <- Sys.time()
```

## System Stats

RAM - 


``` r
get_ram()
```

```
## 90.1 GB
```

CPU - 


``` r
get_cpu()
```

```
## $vendor_id
## [1] "GenuineIntel"
## 
## $model_name
## [1] "Intel(R) Xeon(R) CPU E5-2697 v2 @ 2.70GHz"
## 
## $no_of_cores
## [1] 24
```


Run on Intel(R) Xeon(R) CPU E5-2697 v2 @ 2.70GHz with 24 cores and 9.0146152\times 10^{10} of RAM using Linux version Linux release - 6.10.14-linuxkit and  R version 4.4.0 (2024-04-24)

## Timing

GSEA started at 2025-04-30 14:47:06.577797

GSEA finished at 2025-04-30 14:53:34.127342

GSEA total running time - 


``` r
end_time - start_time
```

```
## Time difference of 6.459159 mins
```

