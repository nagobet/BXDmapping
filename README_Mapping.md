---
title: README_Mapping
author: Nastassia Gobet
date: '2020-10-19'
output:
  html_document:
    keep_md: yes
    code_folding: hide
    fig_caption: yes
    highlight: tango
    number_sections: no
    theme: readable
    toc: yes
    toc_depth: 4
    df_print: paged
  github_document:
    toc: true
    toc_depth: 4
  pdf_document:
    fig_caption: yes
  word_document:
    toc: yes
    toc_depth: '4'
always_allow_html: true
---

# Goal

Highlight important scripts and clarify workflows.




# Scripts

## Coding habits

Keywords used in script names:

* "main": refers to meta-scripts, grouping multiple tasks and calling other scripts.
* "listcommands": refers to list of commands that are then executed (in parallel).
* "OnCluster": refers to scripts adapted to UNIL Wally cluster environnment; can include changes in paths and handling of software stack and SLURM job submission system.

Scripts normally have a line starting by '# GOAL:' to state the aim of the script, except if they are in the "listcommands" category.

Line starting by '##' indicates historical code (not run but may be useful to keep for future).

Line starting by '###' indicates a problem.

## Important scripts

Important (low-level) scripts and what they do.


```r
# create function to retrieve goal from script
retrieveGoal <- function(filename){
  FileInput <- readLines(filename) 
  goal <- grep("# GOAL", FileInput, value=TRUE)
  return(c(filename, goal))
}

# key scripts with goal specified
retrieveGoal("scripts/mapping/MappingOptimization_main.sh")
```

```
## [1] "scripts/mapping/MappingOptimization_main.sh"                                       
## [2] "# GOAL: run mapping various STAR mapping parameters to identify optimal parameters"
```

```r
retrieveGoal("scripts/mapping/MappingOptimizationOnCluster_main.sh")
```

```
## [1] "scripts/mapping/MappingOptimizationOnCluster_main.sh"                              
## [2] "# GOAL: run mapping various STAR mapping parameters to identify optimal parameters"
```

```r
retrieveGoal("scripts/mapping/MergeCount.sh")
```

```
## [1] "scripts/mapping/MergeCount.sh"                                                                     
## [2] "# GOAL: Merge gene expression counts from the different samples (one file per sample) to one file."
```

```r
retrieveGoal("scripts/mapping/filterSamplesCounts.R")
```

```
## [1] "scripts/mapping/filterSamplesCounts.R"                                                              
## [2] "# GOAL: filter out old BXD (05, 29, 29t, 32), problematic BXD (63), parental, and F1 samples counts"
```

```r
retrieveGoal("scripts/mapping/normalizeGeneCountsCPM.R")
```

```
## [1] "scripts/mapping/normalizeGeneCountsCPM.R"
## [2] "# GOAL: normalize gene counts (cpm)"
```

```r
retrieveGoal("scripts/mapping/testingFiltering.R")
```

```
## [1] "scripts/mapping/testingFiltering.R"                                            
## [2] "# GOAL: Perform different filtering of lowly expressed genes and normalization"
```

```r
retrieveGoal("scripts/mapping/eQTLdetectionFilteringTesting_main.sh")
```

```
## [1] "scripts/mapping/eQTLdetectionFilteringTesting_main.sh"               
## [2] "# GOAL: run cis eQTL detection with FastQTL for filtering evaluation"
```

```r
retrieveGoal("scripts/mapping/correctMultiPhenotypeseQTL.R")
```

```
## [1] "scripts/mapping/correctMultiPhenotypeseQTL.R"                                                       
## [2] "# GOAL: compute q-value for eQTL to assess that multiple phenotypes (genes expression) were tested."
```

```r
retrieveGoal("scripts/mapping/BXDspecific_references_databasebased.sh")
```

```
## [1] "scripts/mapping/BXDspecific_references_databasebased.sh"                                   
## [2] "# GOAL: test different personalized references with D2-specific SNPs and indels from dbSNP"
```

```r
retrieveGoal("scripts/mapping/BXDspecific_references_genotypebased.sh")
```

```
## [1] "scripts/mapping/BXDspecific_references_genotypebased.sh"                   
## [2] "# GOAL: test different personalized references using GeneNetwork genotypes"
```

```r
retrieveGoal("scripts/mapping/GenotypeImputation_D2blockmethod.sh")
```

```
## [1] "scripts/mapping/GenotypeImputation_D2blockmethod.sh"              
## [2] "# GOAL: impute genotypes with my D2 blocks homemade simple method"
```

```r
retrieveGoal("scripts/lowlayer/ExtractD2blocksforBXD_main.R")
```

```
## [1] "scripts/lowlayer/ExtractD2blocksforBXD_main.R"
```

```r
retrieveGoal("scripts/lowlayer/extractD2Blocks.R")
```

```
## [1] "scripts/lowlayer/extractD2Blocks.R"
```

```r
# supplementary analyses
retrieveGoal("scripts/mapping/STARmappingB6mm9_main.sh")
```

```
## [1] "scripts/mapping/STARmappingB6mm9_main.sh"                                                                
## [2] "# GOAL: re-run mapping on mm9 non-masked genome+annotation with STAR default parameters (for Figure 2A)."
```

```r
retrieveGoal("scripts/mapping/STARmappingB6mm10primaryassembly_main.sh")
```

```
## [1] "scripts/mapping/STARmappingB6mm10primaryassembly_main.sh"                           
## [2] "# GOAL: re-run mapping on mm10 primary assembly with default parameters and 1-pass."
```

```r
retrieveGoal("scripts/mapping/STARmappingD2_main.sh")
```

```
## [1] "scripts/mapping/STARmappingD2_main.sh"                                              
## [2] "# GOAL: re-run mapping on mm10 primary assembly with default parameters and 1-pass."
```

```r
retrieveGoal("scripts/mapping/STARmappingD2fromB6_main.sh")
```

```
## [1] "scripts/mapping/STARmappingD2fromB6_main.sh"                                                                                                             
## [2] "# GOAL: re-run mapping on D2 reference (customized from B6 mm10 primary assembly, using indels and SNVs from dbSNP)  with default parameters and 1-pass."
```

```r
retrieveGoal("scripts/mapping/genesBiotypes.sh")
```

```
## [1] "scripts/mapping/genesBiotypes.sh"                                                                                              
## [2] "# GOAL: describe the proportion of DM genes in the different genes subtypes. (One of Olivier suggestion in my midthesis exam.)"
## [3] "# GOAL: re-run biotypes analysis with D2 annotation"
```

```r
retrieveGoal("scripts/mapping/plot_biotypes.R")
```

```
## [1] "scripts/mapping/plot_biotypes.R" "# GOAL: plot more biotypes"
```

```r
retrieveGoal("scripts/mapping/mappingKallisto.sh")
```

```
## [1] "scripts/mapping/mappingKallisto.sh"     
## [2] "# GOAL: Mapping with Kallisto in local "
```

```r
retrieveGoal("scripts/mapping/retrieveGenePosition.R")
```

```
## [1] "scripts/mapping/retrieveGenePosition.R"      
## [2] "# GOAL: retrieve gene position from gtf file"
```

```r
retrieveGoal("scripts/mapping/normalizeGeneCountsTPM.R")
```

```
## [1] "scripts/mapping/normalizeGeneCountsTPM.R"
## [2] "# GOAL: normalize gene counts (tpm)"
```

```r
retrieveGoal("scripts/mapping/normalizeGeneCountsTPM_B6mm10.R")
```

```
## [1] "scripts/mapping/normalizeGeneCountsTPM_B6mm10.R"
## [2] "# GOAL: normalize gene counts (tpm)"
```

```r
retrieveGoal("scripts/mapping/periodicity.R")
```

```
## [1] "scripts/mapping/periodicity.R"
```

## Analyses

The higher level analysis and graph production is in analysis/Mapping, with mostly Rmarkdown documents.


## R packages loaded


```r
# retrieve names
myfiles <- c(list.files(path="Mapping", pattern=".[Rmd]$", full.names=TRUE),list.files(path="scripts/mapping", pattern=".[Rmd]$", full.names=TRUE)) 

all_lib <- c()
for(i in myfiles){
  all_lib <- c(all_lib, grep("library\\(", readLines(i), value=TRUE))
}

# clean (remove duplicates)
unique(sort(all_lib))
```

```
## [1] "library(caroline)" "library(edgeR)"    "library(knitr)"   
## [4] "library(limma)"    "library(magick)"   "library(qvalue)"  
## [7] "library(sm)"       "library(vioplot)"  "library(zoo)"
```


# Workflows

The legend used is workflow schemes:


```r
grViz("
digraph legend {

  # a 'graph' statement
  graph [fontsize=16 rankdir=LR]

  # file nodes
  node [shape=box]
  'file (format)'
  
  # script edges
  ' '->'' [label='script']
}
")
```

<!--html_preserve--><div id="htmlwidget-28aa1b1e8d44f073a678" style="width:240px;height:240px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-28aa1b1e8d44f073a678">{"x":{"diagram":"\ndigraph legend {\n\n  # a \"graph\" statement\n  graph [fontsize=16 rankdir=LR]\n\n  # file nodes\n  node [shape=box]\n  \"file (format)\"\n  \n  # script edges\n  \" \"->\"\" [label=\"script\"]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## cis eQTL analysis

```r
grViz("
digraph ciseQTLanalysis {

  # a 'graph' statement
  graph [rankdir=LR fontsize=16]

  # file nodes
  node [shape=box]
  'transcriptome annotation (.gtf)'; 'genotypes (.geno)'; 'normalized gene expression (.tab)';
  'cis eQTL (.txt)'; 'cis eQTL pcorrected (.txt)'

  # script edges
  'transcriptome annotation (.gtf)' -> 'transcriptome annotation (.bed)' [label='create_BED_UCSC.py', fontname=Helvetica]
  'genotypes (.geno)' -> 'genotypes (.vcf)' [label='GeneToVcf.py', fontname=Helvetica]
  'cis eQTL (.txt)' -> 'cis eQTL pcorrected (.txt)' [label='correctMultiPhenotypeseQTL.R', fontname=Helvetica]
  {'genotypes (.vcf)' 'transcriptome annotation (.bed)' 'normalized gene expression (.tab)'} -> 'cis eQTL (.txt)' [label='ciseQTLAnalysis.sh', fontname=Helvetica]
}
")
```

<!--html_preserve--><div id="htmlwidget-4421ecdc55a64322d2f2" style="width:960px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-4421ecdc55a64322d2f2">{"x":{"diagram":"\ndigraph ciseQTLanalysis {\n\n  # a \"graph\" statement\n  graph [rankdir=LR fontsize=16]\n\n  # file nodes\n  node [shape=box]\n  \"transcriptome annotation (.gtf)\"; \"genotypes (.geno)\"; \"normalized gene expression (.tab)\";\n  \"cis eQTL (.txt)\"; \"cis eQTL pcorrected (.txt)\"\n\n  # script edges\n  \"transcriptome annotation (.gtf)\" -> \"transcriptome annotation (.bed)\" [label=\"create_BED_UCSC.py\", fontname=Helvetica]\n  \"genotypes (.geno)\" -> \"genotypes (.vcf)\" [label=\"GeneToVcf.py\", fontname=Helvetica]\n  \"cis eQTL (.txt)\" -> \"cis eQTL pcorrected (.txt)\" [label=\"correctMultiPhenotypeseQTL.R\", fontname=Helvetica]\n  {\"genotypes (.vcf)\" \"transcriptome annotation (.bed)\" \"normalized gene expression (.tab)\"} -> \"cis eQTL (.txt)\" [label=\"ciseQTLAnalysis.sh\", fontname=Helvetica]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## Gene normalization


```r
grViz("
digraph GeneNormalization {

  # a 'graph' statement
  graph [rankdir=LR fontsize=16]

  # file nodes
  node [shape=box]
  'sample counts (.tab)'; 'grouped counts (.tab)'; 'normalized gene expression (.tab)';

  # script edges
  'sample counts (.tab)' -> 'grouped counts (.tab)' [label='MergeCount.sh', fontname=Helvetica]
  'grouped counts (.tab)' -> 'normalized gene expression (.tab)' [label='normalizeGeneCountsCPM.R', fontname=Helvetica]
}
")
```

<!--html_preserve--><div id="htmlwidget-6aae51442c221da11621" style="width:960px;height:288px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-6aae51442c221da11621">{"x":{"diagram":"\ndigraph GeneNormalization {\n\n  # a \"graph\" statement\n  graph [rankdir=LR fontsize=16]\n\n  # file nodes\n  node [shape=box]\n  \"sample counts (.tab)\"; \"grouped counts (.tab)\"; \"normalized gene expression (.tab)\";\n\n  # script edges\n  \"sample counts (.tab)\" -> \"grouped counts (.tab)\" [label=\"MergeCount.sh\", fontname=Helvetica]\n  \"grouped counts (.tab)\" -> \"normalized gene expression (.tab)\" [label=\"normalizeGeneCountsCPM.R\", fontname=Helvetica]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## Filtering evaluation


```r
grViz("
digraph FilteringEvaluation {

  # a 'graph' statement
  graph [rankdir=LR fontsize=16]

  # file nodes
  node [shape=box]
  'grouped counts (.tab)'; 'filtered and normalized gene expression (.tab)'; 'cis eQTL (.txt)'

  # script edges
  'grouped counts (.tab)' -> 'filtered and normalized gene expression (.tab)' [label='testingFiltering.R', fontname=Helvetica]
  'filtered and normalized gene expression (.tab)' -> 'cis eQTL (.txt)' [label='eQTLdetectionFilteringTesting_main.sh', fontname=Helvetica]
}
")
```

<!--html_preserve--><div id="htmlwidget-c6872ae7ec4de72b8019" style="width:960px;height:288px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-c6872ae7ec4de72b8019">{"x":{"diagram":"\ndigraph FilteringEvaluation {\n\n  # a \"graph\" statement\n  graph [rankdir=LR fontsize=16]\n\n  # file nodes\n  node [shape=box]\n  \"grouped counts (.tab)\"; \"filtered and normalized gene expression (.tab)\"; \"cis eQTL (.txt)\"\n\n  # script edges\n  \"grouped counts (.tab)\" -> \"filtered and normalized gene expression (.tab)\" [label=\"testingFiltering.R\", fontname=Helvetica]\n  \"filtered and normalized gene expression (.tab)\" -> \"cis eQTL (.txt)\" [label=\"eQTLdetectionFilteringTesting_main.sh\", fontname=Helvetica]\n}\n","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


# R session information


```r
sessionInfo()
```

```
## R version 3.5.3 (2019-03-11)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=French_Switzerland.1252  LC_CTYPE=French_Switzerland.1252   
## [3] LC_MONETARY=French_Switzerland.1252 LC_NUMERIC=C                       
## [5] LC_TIME=French_Switzerland.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] DiagrammeR_1.0.6.1
## 
## loaded via a namespace (and not attached):
##  [1] visNetwork_2.0.9   digest_0.6.25      jsonlite_1.7.1     magrittr_1.5      
##  [5] evaluate_0.14      rlang_0.4.7        stringi_1.4.6      rstudioapi_0.11   
##  [9] rmarkdown_2.4      RColorBrewer_1.1-2 tools_3.5.3        stringr_1.4.0     
## [13] glue_1.4.2         htmlwidgets_1.5.2  xfun_0.18          yaml_2.2.1        
## [17] compiler_3.5.3     htmltools_0.5.0    knitr_1.30
```

