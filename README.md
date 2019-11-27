HLA Genotype Imputation with Attribute Bagging
======

Kernel Version: v1.4

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

[![Availability](http://www.bioconductor.org/shields/availability/release/HIBAG.svg)](http://www.bioconductor.org/packages/release/bioc/html/HIBAG.html)
[![Years-in-BioC](http://www.bioconductor.org/shields/years-in-bioc/HIBAG.svg)](http://www.bioconductor.org/packages/release/bioc/html/HIBAG.html)
[![Build Status](https://travis-ci.org/zhengxwen/HIBAG.png)](https://travis-ci.org/zhengxwen/HIBAG)
[![Build status](https://ci.appveyor.com/api/projects/status/v650qe8ap4bojxuf?svg=true)](https://ci.appveyor.com/project/zhengxwen/hibag)
[![codecov.io](https://codecov.io/github/zhengxwen/HIBAG/coverage.svg?branch=master)](https://codecov.io/github/zhengxwen/HIBAG?branch=master)


## Features

HIBAG is a state of the art software package for imputing HLA types using SNP data, and it relies on a training set of HLA and SNP genotypes. HIBAG can be used by researchers with published parameter estimates instead of requiring access to large training sample datasets. It combines the concepts of attribute bagging, an ensemble classifier method, with haplotype inference for SNPs and HLA types. Attribute bagging is a technique which improves the accuracy and stability of classifier ensembles using bootstrap aggregating and random variable selection.


## Bioconductor Package

Release Version: 1.22.0

[http://www.bioconductor.org/packages/HIBAG/](http://www.bioconductor.org/packages/HIBAG/)


### Changes in Bioconductor Version (since v1.14.0, Y2017):

* Kernel Version: v1.4
* Modify the kernel to support the GPU extension
* Develop a complementary R package ([HIBAG.gpu](https://github.com/zhengxwen/HIBAG.gpu)) for GPU computing
* The kernel v1.4 outputs exactly the same parameter estimates as v1.3, and the model training with v1.4 is 1.2 times faster than v1.3.


### Changes in Bioconductor Version (since v1.3.0, Y2013):

* Kernel Version: v1.3
* Optimize the calculation of hamming distance using SSE2 and hardware POPCNT instructions if available
* Hardware POPCNT: 2.4x speedup for large-scale data, compared to the implementation in v1.2.4
* SSE2 popcount implementation without hardware POPCNT: 1.5x speedup for large-scale data, compared to the implementation in v1.2.4


## Package Author & Maintainer

Dr. Xiuwen Zheng ([zhengx@u.washington.edu](zhengx@u.washington.edu))


## Pre-fit Model Download

* [https://hibag.s3.amazonaws.com/index.html](https://hibag.s3.amazonaws.com/index.html)
* Platform-specific HLARES models: [https://hibag.s3.amazonaws.com/hlares_index.html](https://hibag.s3.amazonaws.com/hlares_index.html)
* [http://www.biostat.washington.edu/~bsweir/HIBAG/](http://www.biostat.washington.edu/~bsweir/HIBAG/)


## Citation

Zheng, X. *et al*. HIBAG-HLA genotype imputation with attribute bagging. *Pharmacogenomics Journal* 14, 192-200 (2014).
[http://dx.doi.org/10.1038/tpj.2013.18](http://dx.doi.org/10.1038/tpj.2013.18)

Zheng, X. (2018) Imputation-Based HLA Typing with SNPs in GWAS Studies. In: Boegel S. (eds) HLA Typing. Methods in Molecular Biology, Vol 1802. Humana Press, New York, NY. [https://doi.org/10.1007/978-1-4939-8546-3_11](https://doi.org/10.1007/978-1-4939-8546-3_11)


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
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.

* Install the package from the source code:
[download the source code](https://github.com/zhengxwen/HIBAG/tarball/master)
```sh
wget --no-check-certificate https://github.com/zhengxwen/HIBAG/tarball/master -O HIBAG_latest.tar.gz
## or ##
curl -L https://github.com/zhengxwen/HIBAG/tarball/master/ -o HIBAG_latest.tar.gz

## Install ##
R CMD INSTALL HIBAG_latest.tar.gz
```


## Acceleration

### CPU with Intel Intrinsics

* Install the package from the source code with the support of hardware POPCNT (requiring SSE4.2):
You have to customize the package compilation, see: [CRAN: Customizing-package-compilation](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Customizing-package-compilation)

Change `~/.R/Makevars` to, if your machine supports SSE4.2 or higher, assuming GNU Compilers (gcc/g++) or Clang compiler (clang++) are installed:
```sh
## for C code
CFLAGS=-g -O3 -march=native -mtune=native
## for C++ code
CXXFLAGS=-g -O3 -march=native -mtune=native
```
Or force to create hardware POPCNT code:
```sh
## for C code
CFLAGS=-g -O3 -mpopcnt -msse4.2
## for C++ code
CXXFLAGS=-g -O3 -mpopcnt -msse4.2
```

If the package compilation succeeds with hardware POPCNT instructions, you should see a welcome message after loading the package:
```
HIBAG (HLA Genotype Imputation with Attribute Bagging)
Kernel Version: v1.4
Supported by Streaming SIMD Extensions (SSE2 + POPCNT)
```

### GPU with OpenCL

* Install [HIBAG.gpu](https://github.com/zhengxwen/HIBAG.gpu) from Github (for developers/testers only):
```R
library("devtools")
install_github("zhengxwen/HIBAG.gpu")
```
Please use `hlaAttrBagging_gpu()` and `hlaPredict_gpu()` for model training and prediction.


**Speed-up factors for training HIBAG models:**

| CPU (1 core) | CPU (1 core, POPCNT) |
|:------------:|:--------------------:|
| 1            | 1.63 x               |

| 1x NVIDIA Tesla K80 | 1x NVIDIA Tesla M60 | 1x NVIDIA GTX 1080Ti | 1x NVIDIA Tesla P100 | 1x NVIDIA Tesla V100 |
|:-------------------:|:-------------------:|:--------------------:|:--------------------:|:--------------------:|
| 46.5 x              | 57.5 x              | 93.7 x               | 209.1 x              | 246.3 x              |

*using HIBAG v1.14.0 and HIBAG.gpu v0.9.1*

*CPU (1 core), the default installation from Bioconductor supporting SIMD SSE2 instructions, using Intel(R) Xeon(R) CPU E5-2630L @2.40GHz*

*CPU (1 core, POPCNT), optimization with Intel/AMD POPCNT instruction, using Intel(R) Xeon(R) CPU E5-2630L @2.40GHz*

*The benchmark was made possible, in part, through HPC time donated by Microway, Inc. We gratefully acknowledge Microway for providing access to their GPU-accelerated compute cluster (http://www.microway.com/gpu-test-drive/).*


## Archive

[https://github.com/zhengxwen/Archive/tree/master/HIBAG](https://github.com/zhengxwen/Archive/tree/master/HIBAG)
