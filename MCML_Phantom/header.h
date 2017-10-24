#include "cuda_runtime.h"				//#include <cutil.h>
#include "device_launch_parameters.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>

using namespace std;

// DEFINES

#define NUM_BLOCKS 56 //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
//The register usage varies with platform. 64-bit Linux and 32.bit Windows XP have been tested.

#ifdef __linux__ //uses 25 registers per thread (64-bit)
	#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS 17920
#endif

#ifdef _WIN32 //uses 26 registers per thread
	#define NUM_THREADS_PER_BLOCK 288 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
	#define NUM_THREADS 16128
#endif




#define NUMSTEPS_GPU       10000
#define PI                 3.141592654f
#define RPI                0.318309886f
#define MAX_LAYERS         100
#define STR_LEN            200
#define NORMAL             1                    // 1: normal, 0: oblique
#define NUM_DETECTOR    (NORMAL ? 6:9)       // normal: 3 fibers, oblique: 9 fibers
#define ANGLE              (NORMAL ? 0:45)      // normal: 0 degree, oblique: 45 degree
#define SOURCE_NA         (NORMAL ? 0.26:0.26)  // normal: 0.4, oblique; 0.22
#define DETECTOR_NA       (NORMAL ? 0.26:0.26)  // normal: 0.4, oblique; 0.22
#define N_DETECTOR         1.457f
#define N_SOURCE           1.457f
#define ILLUMINATION_R     0.01
#define COLLECT_R          0.01
#define NUMBER_PHOTON      8000000
#define NUMBER_SIMULATION  44
// directory address, maybe needs to be modified
static int wavelength[50] = {400,410,414,420,430,440,450,460,470,480,
490,500,510,520,530,540,542,550,560,570,576,580,590,600,610,620,
630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800}; //can be modified


#define WEIGHTI 0.0001f
//#define WEIGHTI 429497u //0xFFFFFFFFu*WEIGHT
#define CHANCE 0.1f


// TYPEDEFS
typedef struct __align__(16)
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float g;			// Anisotropy factor [-]
	float n;			// Refractive index [-]
}LayerStruct;

typedef struct __align__(16)
{
	float x;				// Global x coordinate [cm]
	float y;				// Global y coordinate [cm]
	float z;				// Global z coordinate [cm]
	float dx;				// (Global, normalized) x-direction
	float dy;				// (Global, normalized) y-direction
	float dz;				// (Global, normalized) z-direction
	float weight;			// Photon weight
	int layer;				// Current layer
	int Ab_z;				// absorbtion is zero
}PhotonStruct;

typedef struct
{
	unsigned long number_of_photons;
	unsigned int n_layers;
	float start_weight;
	LayerStruct* layers;
}SimulationStruct;

typedef struct
{
	float radius[13];
	float NA[13];
	float position[13];
	float angle[13];
	float data[13];
}Fibers;

typedef struct
{
	Fibers* f;
	PhotonStruct* p;						// Pointer to structure array containing all the photon data
	unsigned int* thread_active;			// Pointer to the array containing the thread active status
	unsigned int* num_terminated_photons;	//Pointer to a scalar keeping track of the number of terminated photons
	curandState*  state;
}MemStruct;

typedef struct
{
	float* all;
	float* prob;
	float* cumf;
}G_Array;




__device__ __constant__ unsigned int num_photons_dc[1];
__device__ __constant__ unsigned int n_layers_dc[1];
__device__ __constant__ float start_weight_dc[1];
__device__ __constant__ LayerStruct layers_dc[MAX_LAYERS];


//function
void freeSimulationStruct(SimulationStruct* sim, int n_simulations);
int readData(SimulationStruct** simulations);
void doOneSimulation(SimulationStruct* simulation, int index,ofstream& myfile);

int initMemStruct(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim);
void freeMemStruct(MemStruct* HostMem, MemStruct* DeviceMem);
__global__ void MCd(MemStruct DeviceMem, G_Array DeviceG, unsigned long long seed);
//__global__ void launchPhoton_Global(MemStruct DeviceMem);
int initDCMem(SimulationStruct* sim);
int Write_Simulation_Results(MemStruct* HostMem, SimulationStruct* sim, clock_t simulation_time);
int read_simulation_data(char* filename, SimulationStruct** simulations, int ignoreAdetection);
int interpret_arg(int argc, char* argv[], unsigned long long* seed, int* ignoreAdetection);

__device__ void launchPhoton(PhotonStruct* p, curandState *state);
__global__ void launchPhoton_Global(MemStruct DeviceMem, unsigned long long seed);
__device__ void spin(PhotonStruct* p, G_Array g, curandState *state);
__device__ unsigned int reflect(PhotonStruct*, int, curandState* state);
__device__ unsigned int photonSurvive(PhotonStruct*, curandState* state);
__device__ void atomicAddULL(unsigned long long* address, unsigned int add);
__device__ void detect(PhotonStruct* p, Fibers* f);
__device__ int binarySearch(float *data, float value);
void initFiber(Fibers* f);
void outputFiber(SimulationStruct* sim, float* reflectance, ofstream& myfile);
void calculateReflectance(Fibers* f, float *result);
void inputG(int index, G_Array *g);
int initG(G_Array* HostG, G_Array* DeviceG, int index);
void freeG(G_Array* HostG, G_Array* DeviceG);
