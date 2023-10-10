HLA Genotype Imputation with Attribute Bagging
======

Kernel Version: 1.5

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

[![Availability](http://www.bioconductor.org/shields/availability/release/HIBAG.svg)](http://www.bioconductor.org/packages/release/bioc/html/HIBAG.html)
[![Years-in-BioC](http://www.bioconductor.org/shields/years-in-bioc/HIBAG.svg)](http://www.bioconductor.org/packages/release/bioc/html/HIBAG.html)
[![R](https://github.com/zhengxwen/HIBAG/actions/workflows/r.yml/badge.svg)](https://github.com/zhengxwen/HIBAG/actions/workflows/r.yml)


## Features

HIBAG is a state of the art software package for imputing HLA types using SNP data, and it relies on a training set of HLA and SNP genotypes. HIBAG can be used by researchers with published parameter estimates instead of requiring access to large training sample datasets. It combines the concepts of attribute bagging, an ensemble classifier method, with haplotype inference for SNPs and HLA types. Attribute bagging is a technique which improves the accuracy and stability of classifier ensembles using bootstrap aggregating and random variable selection.


## Bioconductor Package

Release Version: 1.38.0

[http://www.bioconductor.org/packages/HIBAG/](http://www.bioconductor.org/packages/HIBAG/)


### Changes in Bioconductor Version (since v1.26.0, Y2020):

* Kernel Version: **v1.5**
* The kernel v1.5 generates the same training model as v1.4, but 2-6x faster, by taking advantage of Intel AVX, AVX2 and AVX512 intrinsics if available


### Changes in Bioconductor Version (since v1.14.0, Y2017):

* Kernel Version: **v1.4**
* The kernel v1.4 outputs exactly the same model parameter estimates as v1.3, and the model training with v1.4 is 1.2 times faster than v1.3
* Modify the kernel to support the GPU extension


### Changes in Bioconductor Version (since v1.3.0, Y2013):

* Kernel Version: **v1.3**
* Optimize the calculation of hamming distance using SSE2 and hardware POPCNT instructions if available
* Hardware POPCNT: 2.4x speedup for large-scale data, compared to the implementation in v1.2.4
* SSE2 popcount implementation without hardware POPCNT: 1.5x speedup for large-scale data, compared to the implementation in v1.2.4


## Package Author & Maintainer

Dr. Xiuwen Zheng


## Pre-fit Model Download

* [https://hibag.s3.amazonaws.com/index.html](https://hibag.s3.amazonaws.com/index.html)
* Platform-specific HLARES models: [https://hibag.s3.amazonaws.com/hlares_index.html](https://hibag.s3.amazonaws.com/hlares_index.html)


## Citation

Zheng, X. *et al*. HIBAG-HLA genotype imputation with attribute bagging. *Pharmacogenomics Journal* 14, 192-200 (2014).
[doi: 10.1038/tpj.2013.18](http://dx.doi.org/10.1038/tpj.2013.18)

Zheng, X. (2018) Imputation-Based HLA Typing with SNPs in GWAS Studies. In: Boegel S. (eds) HLA Typing. Methods in Molecular Biology, Vol 1802. Humana Press, New York, NY. [doi: 10.1007/978-1-4939-8546-3_11](https://doi.org/10.1007/978-1-4939-8546-3_11)


## Installation

* Bioconductor repository:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("HIBAG")
```

* Development version from Github (for developers/testers only):
```R
library("devtools")
install_github("zhengxwen/HIBAG")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](https://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.


## Acceleration

### CPU with Intel Intrinsics

* **GCC (>= v6.0)** is strongly recommended to compile the HIBAG package (Intel ICC is not suggested).

* `HIBAG::hlaSetKernelTarget("max")` can be used to maximize the algorithm efficiency.


### GPU with OpenCL

* [HIBAG.gpu](https://github.com/zhengxwen/HIBAG.gpu), requiring HIBAG (>= v1.28.0).


## Archive

[https://github.com/zhengxwen/Archive/tree/master/HIBAG](https://github.com/zhengxwen/Archive/tree/master/HIBAG)

[https://bioconductor.org/about/release-announcements](https://bioconductor.org/about/release-announcements)
