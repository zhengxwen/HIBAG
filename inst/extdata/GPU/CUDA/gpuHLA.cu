// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// gpuHLA.cpp: GPU supports for HLA Genotype Imputation
//
// Copyright (C) 2013	Xiuwen Zheng (zhengx@u.washington.edu)
//
// This file is part of HIBAG package.
//
// HIBAG is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// HIBAG is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HIBAG.
// If not, see <http://www.gnu.org/licenses/>.


// g++ -arch x86_64   -I/Developer/NVIDIA/CUDA-5.0/include -I. -I.. -I../../common/inc -o deviceQuery.o -c deviceQuery.cpp
// g++ -arch x86_64  -o deviceQuery deviceQuery.o -Xlinker -rpath /Developer/NVIDIA/CUDA-5.0/lib -L/Developer/NVIDIA/CUDA-5.0/lib -framework CUDA -lcudart 
// mkdir -p ../../bin/darwin/release
// cp deviceQuery ../../bin/darwin/release


// std::system includes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <StructHLA.h>
#include <math.h>

// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

typedef UINT8 TPackedHammingDistance[256u][256u];

__device__ static inline int _HamDist(int Length, UINT8 *SNP,
	UINT8 *H1, UINT8 *H2, TPackedHammingDistance *_PackedHamDist)
{
	int rv = 0;

	for (; Length >= 4; Length -= 4)  // one byte
	{
		rv += (*_PackedHamDist)[*SNP++][(*H1++) | (*H2++ << 1)];
	}
	if (Length > 0)
	{
		UINT8 mask = ~(0xFF << (((UINT8)Length)*2));
		rv += (*_PackedHamDist)[*SNP & mask][((*H1) | (*H2 << 1)) & mask];
	}

	return rv;
}

template<typename TFLOAT, typename TFLOAT_HAPLO>
__global__ void kernal_OutOfBagAcc_F32(int nHLA, int *_HLA_HapIdx, int nHaplo,
	TFLOAT_HAPLO *_HapList, int nGeno, TGPU_Genotype *_GList,
	int nSNP, UINT8 *_RetVal,
	TPackedHammingDistance *_PackedHamDist, TFLOAT *_RareFreq)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < nGeno)
	{
		_GList += idx;

		struct THLAType Guess;
		Guess.Allele1 = Guess.Allele2 = INT_MIN;

		TFLOAT_HAPLO *i1, *i2;
		TFLOAT max=0, prob;

		for (int h1=0; h1 < nHLA; h1++)
		{
			TFLOAT_HAPLO *L1_begin = _HapList + _HLA_HapIdx[h1*2];
			TFLOAT_HAPLO *L1_end   = _HapList + _HLA_HapIdx[h1*2] +
				_HLA_HapIdx[h1*2 + 1];

			// diag value
			prob = 0;
			for (i1=L1_begin; i1 != L1_end; i1++)
			{
				for (i2=i1; i2 != L1_end; i2++)
				{
					prob += ((i1 != i2) ? (2 * i1->Frequency * i2->Frequency) :
							(i1->Frequency * i2->Frequency)) *
						_RareFreq[_HamDist(nSNP, _GList->PackedSNPs,
							i1->PackedHaplo, i2->PackedHaplo, _PackedHamDist)];
				}
			}
			if (max < prob)
			{
				max = prob;
				Guess.Allele1 = Guess.Allele2 = h1;
			}

			// off-diag value
			for (int h2=h1+1; h2 < nHLA; h2++)
			{
				TFLOAT_HAPLO *L2_begin = _HapList + _HLA_HapIdx[h2*2];
				TFLOAT_HAPLO *L2_end   = _HapList + _HLA_HapIdx[h2*2] +
					_HLA_HapIdx[h2*2 + 1];

				prob = 0;
				for (i1=L1_begin; i1 != L1_end; i1++)
				{
					for (i2=L2_begin; i2 != L2_end; i2++)
					{
						prob += (2 * i1->Frequency * i2->Frequency) *
							_RareFreq[_HamDist(nSNP, _GList->PackedSNPs,
							i1->PackedHaplo, i2->PackedHaplo, _PackedHamDist)];
					}
				}
				if (max < prob)
				{
					max = prob;
					Guess.Allele1 = h1; Guess.Allele2 = h2;
				}
			}
		}

		// return 0, 1 or 2

		int T1 = _GList->HLA.Allele1, T2 = _GList->HLA.Allele2;
		int cnt = 0;
		if ((Guess.Allele1==T1) || (Guess.Allele1==T2))
		{
			cnt = 1;
			if (Guess.Allele1==T1) T1 = -1; else T2 = -1;
		}
		if ((Guess.Allele2==T1) || (Guess.Allele2==T2)) cnt ++;

		_RetVal[idx] = cnt;
	}
}


// C
extern "C"
{

static const char *ErrGPUMeg = "GPU fails!";

#define CUDA_MALLOC(ptr, size, cmd)	\
	err = cudaMalloc((void**)&ptr, size); \
	if (err != cudaSuccess) \
	{ \
		fprintf(stderr, "cudaMalloc returned %s\n", cudaGetErrorString(err)); \
		hlaGPU_Finalize(); \
		cmd; \
		throw ErrGPUMeg; \
	}

#define CUDA_MEM_CPY(d, s, size, type, cmd)	\
	err = cudaMemcpy(d, s, size, type); \
	if (err != cudaSuccess) \
	{ \
		fprintf(stderr, "cudaMemcpy returned %s\n", cudaGetErrorString(err)); \
		hlaGPU_Finalize(); \
		cmd; \
		throw ErrGPUMeg; \
	}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

static TPackedHammingDistance *Dev_PackedHammingDistance = NULL;

static double *Dev_RareFreq_Float64 = NULL;
static float  *Dev_RareFreq_Float32 = NULL;

void hlaGPU_Finalize();


/// Export: initialize GPU computing library for HLA imputation
/// throw an error if fails
void hlaGPU_Initialize(uint8_t *hamdist)
{
	cudaError_t err;

	// ***********************************************************************
	// Hamming Distance
	CUDA_MALLOC(Dev_PackedHammingDistance, 256*256, hlaGPU_Finalize());
	CUDA_MEM_CPY(Dev_PackedHammingDistance, hamdist, 256*256,
		cudaMemcpyHostToDevice, hlaGPU_Finalize());


	// ***********************************************************************
	// tables for rare frequencies
	static const double MIN_RARE_FREQ_F64 = 1e-5;
	static const float  MIN_RARE_FREQ_F32 = 1e-5;
	double _RareFreq_Float64[HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2];
	float  _RareFreq_Float32[HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2];

	for (int i=0; i < HIBAG_MAXNUM_SNP_IN_CLASSIFIER*2; i++)
	{
		_RareFreq_Float64[i] = exp (i * log (MIN_RARE_FREQ_F64));
		_RareFreq_Float32[i] = expf(i * logf(MIN_RARE_FREQ_F32));
	}
	_RareFreq_Float64[0] = 1;
	_RareFreq_Float32[0] = 1;

	CUDA_MALLOC(Dev_RareFreq_Float64, sizeof(_RareFreq_Float64),
		hlaGPU_Finalize());
	CUDA_MEM_CPY(Dev_RareFreq_Float64, _RareFreq_Float64,
		sizeof(_RareFreq_Float64), cudaMemcpyHostToDevice, hlaGPU_Finalize());
	CUDA_MALLOC(Dev_RareFreq_Float32, sizeof(_RareFreq_Float32),
		hlaGPU_Finalize());
	CUDA_MEM_CPY(Dev_RareFreq_Float32, _RareFreq_Float32,
		sizeof(_RareFreq_Float32), cudaMemcpyHostToDevice, hlaGPU_Finalize());
}



/// Export: initialize GPU computing library for HLA imputation
void hlaGPU_Finalize()
{
	if (Dev_PackedHammingDistance)
	{
		cudaFree(Dev_PackedHammingDistance);
		Dev_PackedHammingDistance = NULL;
	}
	if (Dev_RareFreq_Float64)
	{
		cudaFree(Dev_RareFreq_Float64);
		Dev_RareFreq_Float64 = NULL;
	}
	if (Dev_RareFreq_Float32)
	{
		cudaFree(Dev_RareFreq_Float32);
		Dev_RareFreq_Float32 = NULL;
	}
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int hlaGPU_OutOfBagAcc_F32(int nHLA, int *_HLA_HapIdx, int nHaplo,
	TGPU_Haplotype_F32 *_HapList, int nGeno, TGPU_Genotype *_GList, int nSNP)
{
	if (nGeno <= 0) return 0;

	// copy host memory to device memory
	cudaError_t err;

	int *Dev_HLA_HapIdx = NULL;
	TGPU_Haplotype_F32 *Dev_HapList = NULL;
	TGPU_Genotype *Dev_GList = NULL;
	uint8_t *Dev_RetVal = NULL;

	// allocate memory
	CUDA_MALLOC(Dev_HLA_HapIdx, sizeof(int)*nHLA*2, {});
	CUDA_MALLOC(Dev_HapList, sizeof(TGPU_Haplotype_F32)*nHaplo,
		cudaFree(Dev_HLA_HapIdx));
	CUDA_MALLOC(Dev_GList, sizeof(TGPU_Genotype)*nGeno,
		{ cudaFree(Dev_HLA_HapIdx); cudaFree(Dev_HapList); });
	CUDA_MALLOC(Dev_RetVal, sizeof(uint8_t)*nGeno,
		{ cudaFree(Dev_HLA_HapIdx); cudaFree(Dev_HapList); cudaFree(Dev_GList); });

	// memory copy
	CUDA_MEM_CPY(Dev_HLA_HapIdx, _HLA_HapIdx, sizeof(int)*nHLA*2,
		cudaMemcpyHostToDevice,  { cudaFree(Dev_HLA_HapIdx);
		cudaFree(Dev_HapList); cudaFree(Dev_GList); cudaFree(Dev_RetVal); });
	CUDA_MEM_CPY(Dev_HapList, _HapList, sizeof(TGPU_Haplotype_F32)*nHaplo,
		cudaMemcpyHostToDevice,  { cudaFree(Dev_HLA_HapIdx);
		cudaFree(Dev_HapList); cudaFree(Dev_GList); cudaFree(Dev_RetVal); });
	CUDA_MEM_CPY(Dev_GList, _GList, sizeof(TGPU_Genotype)*nGeno,
		cudaMemcpyHostToDevice,  { cudaFree(Dev_HLA_HapIdx);
		cudaFree(Dev_HapList); cudaFree(Dev_GList); cudaFree(Dev_RetVal); });

	// run in parallel
	dim3 dimBlock(32);
	dim3 dimGrid((nGeno/32) + ((nGeno % 32) ? 1:0));

	kernal_OutOfBagAcc_F32<float, TGPU_Haplotype_F32><<<dimGrid, dimBlock>>>(
		nHLA, Dev_HLA_HapIdx, nHaplo, Dev_HapList, nGeno, Dev_GList, nSNP,
		Dev_RetVal, Dev_PackedHammingDistance, Dev_RareFreq_Float32);


	// merge results
	int rv = 0;
	if (nGeno <= 16384)
	{
		UINT8 buffer[16384];
		CUDA_MEM_CPY(buffer, Dev_RetVal, nGeno,
			cudaMemcpyDeviceToHost,  { cudaFree(Dev_HLA_HapIdx);
			cudaFree(Dev_HapList); cudaFree(Dev_GList); cudaFree(Dev_RetVal); });
		for (int i=0; i < nGeno; i++) rv += buffer[i];
	} else {
		UINT8 *buffer = new UINT8[nGeno];
		CUDA_MEM_CPY(buffer, Dev_RetVal, nGeno,
			cudaMemcpyDeviceToHost,  { cudaFree(Dev_HLA_HapIdx);
			cudaFree(Dev_HapList); cudaFree(Dev_GList); cudaFree(Dev_RetVal);
			delete []buffer; });
		for (int i=0; i < nGeno; i++) rv += buffer[i];
		delete []buffer;
	}

	cudaFree(Dev_HLA_HapIdx);
	cudaFree(Dev_HapList);
	cudaFree(Dev_GList);
	cudaFree(Dev_RetVal);

	return rv;
}


int hlaGPU_OutOfBagAcc_F64(int nHLA, int *_HLA_HapIdx,
	TGPU_Haplotype_F64 *_HapList, int nGeno, TGPU_Genotype *_GList)
{
	return 0;	
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#if CUDART_VERSION < 5000

/// This function wraps the CUDA Driver API into a template function
static inline void getCudaAttribute(int *attribute,
	CUdevice_attribute device_attribute, int device)
{
	CUresult error = cuDeviceGetAttribute(attribute, device_attribute, device);
	if (CUDA_SUCCESS != error)
	{
		fprintf(stderr,
			"cuSafeCallNoSync() Driver API error = %04d from file <%s>, line %i.\n",
			error, __FILE__, __LINE__);
		throw ErrGPUMeg;
	}
}

#endif


/// Export: enumerate the properties of the CUDA devices present in the system
/// throw an error if fails
void hlaGPU_DeviceQuery()
{
	printf("CUDA Device Query:\n\n");

	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess)
	{
		fprintf(stderr, "cudaGetDeviceCount returned %d\n-> %s\n",
			(int)error_id, cudaGetErrorString(error_id));
		throw ErrGPUMeg;
	}

	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount <= 0)
	{
		fprintf(stderr, "There are no available device(s) that support CUDA\n");
		throw ErrGPUMeg;
	} else {
		printf("Detected %d CUDA Capable device(s)\n", deviceCount);
	}

	int dev, driverVersion = 0, runtimeVersion = 0;
	for (dev = 0; dev < deviceCount; dev++)
    {
		cudaSetDevice(dev);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);

		printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

		// Console log
		cudaDriverGetVersion(&driverVersion);
		cudaRuntimeGetVersion(&runtimeVersion);
		printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
			driverVersion/1000, (driverVersion%100)/10,
			runtimeVersion/1000, (runtimeVersion%100)/10);
		printf("  CUDA Capability Major/Minor version number:    %d.%d\n",
			deviceProp.major, deviceProp.minor);

		printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
			(float)deviceProp.totalGlobalMem/1048576.0f,
			(unsigned long long) deviceProp.totalGlobalMem);

		printf("  (%2d) Multiprocessors x (%3d) CUDA Cores/MP:    %d CUDA Cores\n",
			deviceProp.multiProcessorCount,
			_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
			_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
			deviceProp.multiProcessorCount);
		printf("  GPU Clock rate:                                %.0f MHz (%0.2f GHz)\n",
			deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);

#if CUDART_VERSION >= 5000
		// This is supported in CUDA 5.0 (runtime API device properties)
		printf("  Memory Clock rate:                             %.0f Mhz\n",
			deviceProp.memoryClockRate * 1e-3f);
		printf("  Memory Bus Width:                              %d-bit\n",
			deviceProp.memoryBusWidth);

		if (deviceProp.l2CacheSize)
		{
			printf("  L2 Cache Size:                                 %d bytes\n",
				deviceProp.l2CacheSize);
		}
#else
		// This only available in CUDA 4.0-4.2 (but these were only exposed in the CUDA Driver API)
		int memoryClock;
		getCudaAttribute(&memoryClock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE, dev);
		printf("  Memory Clock rate:                             %.0f Mhz\n",
			memoryClock * 1e-3f);
		int memBusWidth;
		getCudaAttribute(&memBusWidth, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, dev);
		printf("  Memory Bus Width:                              %d-bit\n",
			memBusWidth);
		int L2CacheSize;
		getCudaAttribute(&L2CacheSize, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, dev);

		if (L2CacheSize)
		{
			printf("  L2 Cache Size:                                 %d bytes\n",
				L2CacheSize);
		}
#endif

		printf("  Max Texture Dimension Size (x,y,z)             1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
			deviceProp.maxTexture1D   , deviceProp.maxTexture2D[0], deviceProp.maxTexture2D[1],
			deviceProp.maxTexture3D[0], deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
        printf("  Max Layered Texture Size (dim) x layers        1D=(%d) x %d, 2D=(%d,%d) x %d\n",
			deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1],
			deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1],
			deviceProp.maxTexture2DLayered[2]);

		printf("  Total amount of constant memory:               %lu bytes\n",
			deviceProp.totalConstMem);
		printf("  Total amount of shared memory per block:       %lu bytes\n",
			deviceProp.sharedMemPerBlock);
		printf("  Total number of registers available per block: %d\n",
			deviceProp.regsPerBlock);
		printf("  Warp size:                                     %d\n",
			deviceProp.warpSize);
		printf("  Maximum number of threads per multiprocessor:  %d\n",
			deviceProp.maxThreadsPerMultiProcessor);
		printf("  Maximum number of threads per block:           %d\n",
			deviceProp.maxThreadsPerBlock);
		printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
			deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
			deviceProp.maxThreadsDim[2]);
		printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
			deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
			deviceProp.maxGridSize[2]);
		printf("  Maximum memory pitch:                          %lu bytes\n",
			deviceProp.memPitch);
		printf("  Texture alignment:                             %lu bytes\n",
			deviceProp.textureAlignment);
		printf("  Concurrent copy and kernel execution:          %s with %d copy engine(s)\n",
			(deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
		printf("  Run time limit on kernels:                     %s\n",
			deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
		printf("  Integrated GPU sharing Host Memory:            %s\n",
			deviceProp.integrated ? "Yes" : "No");
		printf("  Support host page-locked memory mapping:       %s\n",
			deviceProp.canMapHostMemory ? "Yes" : "No");
		printf("  Alignment requirement for Surfaces:            %s\n",
			deviceProp.surfaceAlignment ? "Yes" : "No");
		printf("  Device has ECC support:                        %s\n",
			deviceProp.ECCEnabled ? "Enabled" : "Disabled");
#ifdef WIN32
		printf("  CUDA Device Driver Mode (TCC or WDDM):         %s\n",
			deviceProp.tccDriver ?
			"TCC (Tesla Compute Cluster Driver)" : "WDDM (Windows Display Driver Model)");
#endif
		printf("  Device supports Unified Addressing (UVA):      %s\n",
			deviceProp.unifiedAddressing ? "Yes" : "No");
		printf("  Device PCI Bus ID / PCI location ID:           %d / %d\n",
			deviceProp.pciBusID, deviceProp.pciDeviceID);

		const static char *sComputeMode[] =
		{
			"Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
			"Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
			"Prohibited (no host thread can use ::cudaSetDevice() with this device)",
			"Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
			"Unknown",
			NULL
		};
		printf("  Compute Mode:\n");
		printf("     < %s >\n", sComputeMode[deviceProp.computeMode]);
	}

	// csv masterlog info
	// *****************************
	// exe and CUDA driver name
	printf("\n");
	printf("deviceQuery, CUDA Driver = CUDART");

	// driver version
	printf(", CUDA Driver Version = %d.%d",
		driverVersion/1000, (driverVersion%100)/10);
	// Runtime version
	printf(", CUDA Runtime Version = %d.%d",
		runtimeVersion/1000, (runtimeVersion%100)/10);
	// Device count
	printf(", NumDevs = %d", deviceCount);

	// Print Out all device Names
	for (dev = 0; dev < deviceCount; ++dev)
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);
		printf(", Device%d = %s", dev, deviceProp.name);
	}

	printf("\n");
}


} // extern "C"
