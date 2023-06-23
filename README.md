# Brobustedger
Aim of this package is used to identify differentially expressed genes from RNA-seq count data through **boosting robust edgeR**  when the dataset contains outliers or missing values. Missing values are imputed using random forest method. Finally DEGs are calculated by robust edgeR. 
## Installation
```r
#Install edgeR from bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
# Install devtools from CRAN
install.packages("devtools")

# development version from GitHub:
require("devtools")
devtools::install_github("BandhanSarker/bRobustEdgeR")
```
## Usage

```r
#Identify differentially expressed gene list
DEGList (data, n1, n2, p.threshold = 0.05)
```
## Arguments
* **data**	  RNA-seq count data matrix, whose row contains genes and column contains  samples                 
* **n1**	  control group
* **n2**	  case group
* **p.threshold**	adjusted p-value by default=0.05
## Value 
Differentially expressed genes list
## Examples
```r
 data (sampleDataOut)
 DEGList(sampleDataOut,3,3,0.05)
```

```r
#Sample RNA-seq count data with outliers value
data(sampleDataOut)
```
```r
#Sample RNA-seq count data with missing value
data(sampleDataMiss)
```
 
## Author(s): 
Bandhan Sarker, Md. Matiur Rahaman, Muhammad Habibulla Alamin, Md. Ariful Islam and Md.  Nurul Haque Mollah
