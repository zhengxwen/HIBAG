\name{hlaSetKernelTarget}
\alias{hlaSetKernelTarget}
\title{Set the CPU target}
\description{
    Set the CPU target that the HIBAG algorithm is built on.
}
\usage{
hlaSetKernelTarget(
    cpu=c("auto", "avx2", "avx", "sse4", "sse2", "base", "avx512bw"))
}
\arguments{
    \item{cpu}{Specify the Intel/AMD CPU flag; "auto" by default}
}
\details{
    If \code{cpu="auto"}, the kernel target will be automatically determined
according to the CPU capabilities. By default, "avx2" is used instead of
"avx512bw" even if the CPU supports the AVX512F+AVX512BW intrinsics, since
the CPU may reduce the frequency of the cores dynamically to keep power usage
of AVX512 within bounds.

    The HIBAG algorithm is optimized using SIMD intrinsics on different CPU
platforms to leverage the efficiency of the target Intel/AMD platform. The
higher version of the C++ compiler is needed to enable the compilation of
AVX2 and AVX512F intrinsics, e.g., GCC >= v6.0.
If the compiler does not support the CPU target, the implementation on that
target will be disabled.
}
\value{
    Return a string describing the CPU capabilities.
}
\author{Xiuwen Zheng}

\examples{
hlaSetKernelTarget("auto")
}

\keyword{HLA}
\keyword{CPU}
