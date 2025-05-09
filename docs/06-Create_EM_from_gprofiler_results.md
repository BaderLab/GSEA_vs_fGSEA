# Create Enrichment map from R with g:Profiler results

## Initialize variables and libraries


``` r
#use library
#make sure biocManager is installed
tryCatch(expr = { library("BiocManager")}, 
         error = function(e) { 
           install.packages("BiocManager")}, 
         finally = library("BiocManager"))

tryCatch(expr = { library("ggplot2")}, 
         error = function(e) { install.packages("ggplot2")}, 
         finally = library("ggplot2"))

#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library("RCy3")}, 
         error = function(e) { BiocManager::install("RCy3")}, 
         finally = library("RCy3"))

tryCatch(expr = { library("httr")}, 
         error = function(e) { BiocManager::install("httr")}, 
         finally = library("httr"))
```

## Configurable Parameters


``` r
# is_docker - true/false depending on if you are running R from docker
is_docker <- TRUE

#directory where all the original input data file are
# for example ./data/
working_dir <- "./data/"


#directory where all the generated data files are found.
# For example - ./generated_data/
# If you are using all the notebooks from this set the generated data will be
# put in the ./generated_data folder.  You have to specify if it is gsea or 
# gprofiler
output_dir <- "./generated_data/g_profiler"


#defined threshold for gprofiler enrichments 
#p-value to filter all the genesets.  For example -   1.0
pvalue_gprofiler_threshold <- 1.0

#q-value to filter all the genesets.  For example -   0.05
qvalue_gprofiler_threshold <- 0.001

#similarity threshold to filter all the genesets connections/edges.  
# For example -   0.375
similarity_threshold <- "0.35"

#similarity metric to filter all the genesets connections/edges 
# (can be OVERLAP, JACCARD, or COMBINED.   For example -   Combined
similarity_metric = "JACCARD"
```

## Specify Data files
Depending on whether you are creating your enrichment map from g:Profiler or GSEA results the sets of files might be a little different.  Minimally, you will need to specify:
  * gmt file
  * enrichment results file

We have multiple g:profiler results.  
  * varied geneset size limit ( 250, 1000 or 10,000)
  * varied geneset sources - baderlab genesets or g:profiler sets.


``` r
#get the newest gprofiler gmt file
gprof_files <- file.info(list.files(file.path(getwd(),output_dir),
                                    pattern = "gprofiler_full",
                                    full.names = TRUE) )
gmt_gprofiler_file <- rownames(gprof_files)[which.max(gprof_files$mtime)]
#gmt_gprofiler_file<-file.path(getwd(),output_dir,
#                "gprofiler_full_hsapiens_e109_eg56_p17_1d3191d_.name.gmt")

gprofiler_results_filename1 <-file.path(getwd(),output_dir,
                "gProfiler_hsapiens_lab2_results_GEM_termmin3_max250.gem.txt")

gprofiler_results_filename2 <-file.path(getwd(),output_dir,
                "gProfiler_hsapiens_lab2_results_GEM_termmin3_max10000.gem.txt")


current_network_name <- paste("gprofiler_max250",pvalue_gprofiler_threshold,
                              qvalue_gprofiler_threshold,sep="_")
```


## Launch Cytoscape

Launch Cytoscape (by default cytoscape will automatically enable rest so as long as cytoscape 3.3 or higher is open R should be able to communicate with it).  Make sure if you get an message asking you if you want communicate with other apps that you select "Allow".  

## Make sure you can connect to Cytoscape

``` r
if(is_docker){
  current_base = "host.docker.internal:1234/v1"
  .defaultBaseUrl <- "http://host.docker.internal:1234/v1"
} else{
  current_base = "localhost:1234/v1"
}

cytoscapePing (base.url = current_base)
```

```
## You are connected to Cytoscape!
```

``` r
cytoscapeVersionInfo (base.url = current_base)
```

```
##       apiVersion cytoscapeVersion 
##             "v1"         "3.10.2"
```
***
## Create an Enrichment map

If you are running R from within a docker you need to first upload your datafiles to Cytoscape before you can create your enrichment map


``` r
#if using docker we need to replace all the the paths to the host path
if(is_docker) {
  upload_em_file <- function(localPath) {
    bname <- basename(localPath)
    r <- POST(
      url = 
paste('http://host.docker.internal:1234/enrichmentmap/textfileupload?fileName=', 
                  bname, sep=""),
      config = list(),
      body = list(file = upload_file(localPath)),
      encode = "multipart",
      handle = NULL
    )
    content(r,"parsed")$path
  }
  
  # "upload" the files to the host machine and replace each path with the 
  # host machine path
   gmt_gprofiler_file <- upload_em_file(gmt_gprofiler_file)
  gprofiler_results_filename1 <- upload_em_file(gprofiler_results_filename1)
  gprofiler_results_filename2 <- upload_em_file(gprofiler_results_filename2)
}
```
***
## Create an Enrichment map - run EM command

``` r
#######################################
#create EM

em_command = paste('enrichmentmap build analysisType="generic" gmtFile=',
                   gmt_gprofiler_file,
                   'pvalue=',pvalue_gprofiler_threshold, 
                   'qvalue=',qvalue_gprofiler_threshold,
                   'similaritycutoff=',similarity_threshold,
                   'coefficients=',similarity_metric,
                   'enrichmentsDataset1=',gprofiler_results_filename1, 
                   'enrichmentsDataset2=',gprofiler_results_filename2,
                   'gmtFile=',gmt_gprofiler_file,
                   'filterByExpressions=false',
                   sep=" ")

#enrichment map command will return the suid of newly created network.
response <- commandsGET(em_command,base.url = current_base)

current_network_suid <- 0
#enrichment map command will return the suid of newly created network 
# unless it Failed.  If it failed it will contain the word failed
if(grepl(pattern="Failed", response)){
  paste(response)
} else {
  current_network_suid <- response
}

#check to see if the network name is unique
current_names <- getNetworkList(base.url = current_base)
if(current_network_name %in% current_names){
  #if the name already exists in the network names then put the SUID in front
  # of the name (this does not work if you put the suid at the end of the name)
  current_network_name <- paste(current_network_suid,
                                current_network_name,  sep="_")
}
response <- renameNetwork(title=current_network_name, 
                       network = as.numeric(current_network_suid),
                       base.url = current_base)
```

## Get a screen shot of the initial network.

``` r
#you can only output the file if it isn't on docker
#on docker is put it into the user's home directory with docker has not access to
if(!is_docker){
  output_network_file <- file.path(getwd(),"initial_screenshot_network.png")
  output_network_file_current <- output_network_file

  fitContent()

  if(file.exists(output_network_file)){
    #cytoscape hangs waiting for user response if file already exists.  
    # Remove it first
    response <- file.remove(output_network_file)
  } 

  response <- exportImage(output_network_file, type = "png",
                          base.url = current_base)
}
```


Change the files and create all the different networks we generated in class in Cytoscape.  

