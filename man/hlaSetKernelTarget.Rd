\name{hlaSetKernelTarget}
\alias{hlaSetKernelTarget}
\title{Set the CPU target}
\description{
    Set the CPU target that the HIBAG algorithm is built on.
}
\usage{
hlaSetKernelTarget(cpu=c("max", "auto.avx2",
    "base", "sse2", "sse4", "avx", "avx2", "avx512f", "avx512bw"))
}
\arguments{
    \item{cpu}{Specify the Intel/AMD CPU flag; "max" by default}
}
\details{
    If \code{cpu="max"}, the kernel target will be automatically determined
according to the CPU capabilities to maximize the algorithm efficiency.
When \code{cpu="auto.avx2"}, "avx2" is used instead of "avx512f", "avx512bw"
even if the CPU supports the AVX512F and AVX512BW intrinsics, since the CPU
may reduce the frequency of the cores dynamically to keep power usage of
AVX512 within bounds; if AVX2 is not applicable, other target will be
automatically determined.

    The HIBAG algorithm is optimized using different SIMD instruction sets to
leverage the efficiency of the target Intel/AMD platform. The higher version of
the C++ compiler is needed to enable the compilation of AVX2 and AVX512F
intrinsics, e.g., GCC >= v6.0. If the compiler does not support the CPU target,
the implementation on that target will be disabled.
}
\value{
    Return a character vector for describing the CPU capabilities, the compiler
information and the supported implementation.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{hlaAttrBagging}}, \code{\link{hlaParallelAttrBagging}},
    \code{\link{predict.hlaAttrBagClass}}, \code{\link{hlaPredict}}
}

\examples{
hlaSetKernelTarget("auto")
}

\keyword{HLA}
\keyword{CPU}
