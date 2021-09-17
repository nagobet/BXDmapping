---
title: "Coverage"
author: "Nastassia Gobet"
date: "31 July 2019 -"
output:
  html_document:
    code_folding: hide
    fig_caption: yes
    highlight: tango
    number_sections: no
    theme: readable
    toc: yes
    toc_depth: 4
    toc_float:
      collapsed: true
      smooth_scroll: false
    df_print: paged
    keep_md: true
  pdf_document:
    fig_caption: yes
  word_document:
    toc: yes
    toc_depth: '4'
---




```r
# store all paths to plots
paths <- array(data="", dim=c(8,20))
rownames(paths) <- c("B61nsd", "B62nsd", "DB1nsd", "DB2nsd", "LB61nsd", "LB62nsd", "LDB1nsd", "LDB2nsd")
colnames(paths) <- c(1:19,"X")
##paths

for(chr in colnames(paths)){
  for(s in rownames(paths)){
  paths[s, chr]  <- paste("plots_log_shiftEnd/", s, "_", chr, ".png", sep="")
  }
}
```


# plot by chromosome

The order of samples is B61nsd, B62nsd, DB1nsd, DB2nsd, LB61nsd, LB62nsd, LDB1nsd, LDB2nsd. Biological replicates are on the sample row.


```r
include_graphics(paths[,1])
```

<img src="plots_log_shiftEnd/B61nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_1.png" width="50%" />

```r
include_graphics(paths[,2])
```

<img src="plots_log_shiftEnd/B61nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_2.png" width="50%" />

```r
include_graphics(paths[,3])
```

<img src="plots_log_shiftEnd/B61nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_3.png" width="50%" />

```r
include_graphics(paths[,4])
```

<img src="plots_log_shiftEnd/B61nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_4.png" width="50%" />

```r
include_graphics(paths[,5])
```

<img src="plots_log_shiftEnd/B61nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_5.png" width="50%" />

```r
include_graphics(paths[,6])
```

<img src="plots_log_shiftEnd/B61nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_6.png" width="50%" />

```r
include_graphics(paths[,7])
```

<img src="plots_log_shiftEnd/B61nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_7.png" width="50%" />

```r
include_graphics(paths[,8])
```

<img src="plots_log_shiftEnd/B61nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_8.png" width="50%" />

```r
include_graphics(paths[,9])
```

<img src="plots_log_shiftEnd/B61nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_9.png" width="50%" />

```r
include_graphics(paths[,10])
```

<img src="plots_log_shiftEnd/B61nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_10.png" width="50%" />

```r
include_graphics(paths[,11])
```

<img src="plots_log_shiftEnd/B61nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_11.png" width="50%" />

```r
include_graphics(paths[,12])
```

<img src="plots_log_shiftEnd/B61nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_12.png" width="50%" />

```r
include_graphics(paths[,13])
```

<img src="plots_log_shiftEnd/B61nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_13.png" width="50%" />

```r
include_graphics(paths[,14])
```

<img src="plots_log_shiftEnd/B61nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_14.png" width="50%" />

```r
include_graphics(paths[,15])
```

<img src="plots_log_shiftEnd/B61nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_15.png" width="50%" />

```r
include_graphics(paths[,16])
```

<img src="plots_log_shiftEnd/B61nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_16.png" width="50%" />

```r
include_graphics(paths[,17])
```

<img src="plots_log_shiftEnd/B61nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_17.png" width="50%" />

```r
include_graphics(paths[,18])
```

<img src="plots_log_shiftEnd/B61nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_18.png" width="50%" />

```r
include_graphics(paths[,19])
```

<img src="plots_log_shiftEnd/B61nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_19.png" width="50%" />

```r
include_graphics(paths[,"X"])
```

<img src="plots_log_shiftEnd/B61nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_X.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_X.png" width="50%" />

```r
##Graphs <- c("plots/B61nsd_1.png", "plots/B62nsd_1.png")
##include_graphics(unlist(Graphs))

##include_graphics(c(paths["LB62nsd","X"],paths["LB62nsd","X"],paths["LB62nsd","X"],paths["B62nsd",4],paths["LB62nsd","X"],paths["LB62nsd","X"],paths["LB62nsd","X"],paths["LB62nsd","X"]))
```


# plot by sample


```r
include_graphics(paths["B61nsd",])
```

<img src="plots_log_shiftEnd/B61nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/B61nsd_X.png" width="50%" />

```r
include_graphics(paths["B62nsd",])
```

<img src="plots_log_shiftEnd/B62nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/B62nsd_X.png" width="50%" />

```r
include_graphics(paths["DB1nsd",])
```

<img src="plots_log_shiftEnd/DB1nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/DB1nsd_X.png" width="50%" />

```r
include_graphics(paths["DB2nsd",])
```

<img src="plots_log_shiftEnd/DB2nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/DB2nsd_X.png" width="50%" />

```r
include_graphics(paths["LB61nsd",])
```

<img src="plots_log_shiftEnd/LB61nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LB61nsd_X.png" width="50%" />

```r
include_graphics(paths["LB62nsd",])
```

<img src="plots_log_shiftEnd/LB62nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LB62nsd_X.png" width="50%" />

```r
include_graphics(paths["LDB1nsd",])
```

<img src="plots_log_shiftEnd/LDB1nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LDB1nsd_X.png" width="50%" />

```r
include_graphics(paths["LDB2nsd",])
```

<img src="plots_log_shiftEnd/LDB2nsd_1.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_2.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_3.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_4.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_5.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_6.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_7.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_8.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_9.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_10.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_11.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_12.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_13.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_14.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_15.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_16.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_17.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_18.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_19.png" width="50%" /><img src="plots_log_shiftEnd/LDB2nsd_X.png" width="50%" />


# Session information

This document is built with R version 3.5.1. More info:


```r
# Rstudio session info
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
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
## [1] knitr_1.23
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     htmltools_0.3.6
##  [5] yaml_2.2.0      Rcpp_1.0.2      stringi_1.4.3   rmarkdown_1.14 
##  [9] stringr_1.4.0   xfun_0.8        digest_0.6.20   evaluate_0.14
```
