#include "header.h"
#include <float.h> //for FLT_MAX



__device__ float rnGen(curandState *s)
{
	float x = curand_uniform(s);
    return x;
}

void doOneSimulation(SimulationStruct* simulation, int index, ofstream& myfile)
{
	extern int wavelength[];
	unsigned long long seed = time(NULL);
	float reflectance[8] = {0};

	MemStruct DeviceMem;
	MemStruct HostMem;
	unsigned int threads_active_total=1;
	unsigned int i,ii;

    cudaError_t cudastat;

	initMemStruct(&HostMem,&DeviceMem,simulation);
	initDCMem(simulation);

    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid(NUM_BLOCKS);

	launchPhoton_Global<<<dimGrid,dimBlock>>>(DeviceMem, seed);
	cudaThreadSynchronize();     //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
	cudastat=cudaGetLastError(); // Check if there was an error
	if(cudastat)printf("Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));

	i=0;

	// initialize the g_factor
	G_Array HostG, DeviceG;
	initG(&HostG,&DeviceG,index);

	while(threads_active_total>0)
	{
		i++;
		initFiber(HostMem.f);
	    cudaMemcpy(DeviceMem.f,HostMem.f,NUM_THREADS*sizeof(Fibers),cudaMemcpyHostToDevice);

		//run the kernel
		seed = time(NULL);
		MCd<<<dimGrid,dimBlock>>>(DeviceMem,DeviceG, seed);
		cudaThreadSynchronize(); //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
		cudastat=cudaGetLastError(); // Check if there was an error
		if(cudastat)printf("Error code=%i, %s.\n",cudastat,cudaGetErrorString(cudastat));

		// Copy thread_active from device to host, later deleted
		cudaMemcpy(HostMem.thread_active, DeviceMem.thread_active, NUM_THREADS*sizeof(unsigned int), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.thread_active,DeviceMem.thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyDeviceToHost) );
		threads_active_total = 0;
		for(ii=0;ii<NUM_THREADS;ii++) threads_active_total+=HostMem.thread_active[ii];

		//CUDA_SAFE_CALL(cudaMemcpy(HostMem.num_terminated_photons,DeviceMem.num_terminated_photons,sizeof(unsigned int),cudaMemcpyDeviceToHost) );

		//printf("Run %u, Number of photons terminated %u, Threads active %u\n",i,*HostMem.num_terminated_photons,threads_active_total);

		cudaMemcpy(HostMem.f, DeviceMem.f, NUM_THREADS*sizeof(Fibers), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.f,DeviceMem.f,NUM_THREADS*sizeof(Fibers),cudaMemcpyDeviceToHost));
		calculateReflectance(HostMem.f,reflectance);
	}
	//cout << "#" << index << " Simulation done!\n";

	myfile << wavelength[index] << "\t";
	outputFiber(simulation, reflectance, myfile);
	freeMemStruct(&HostMem,&DeviceMem);
}

void calculateReflectance(Fibers* f, float *result)
{
	for(int i = 0; i < NUM_THREADS; i++)
	{
		if(NORMAL)
		{
			result[0] += f[i].data[1];
			result[1] += f[i].data[2];
			result[2] += f[i].data[3];
			result[3] += f[i].data[4];
			result[4] += f[i].data[5];
			result[5] += f[i].data[6];
		}
		else
		{
			result[0] += f[i].data[1];
			result[1] += f[i].data[2];
			result[2] += f[i].data[3];
			result[3] += f[i].data[4] + f[i].data[7];
			result[4] += f[i].data[5] + f[i].data[8];
			result[5] += f[i].data[6] + f[i].data[9];
		}
	}
}

//Device function to add an unsigned integer to an unsigned long long using CUDA Compute Capability 1.1
__device__ void atomicAddULL(unsigned long long* address, unsigned int add)
{
	if(atomicAdd((unsigned int*)address,add)+add<add)
		atomicAdd(((unsigned int*)address)+1,1u);
}

__global__ void MCd(MemStruct DeviceMem, G_Array DeviceG, unsigned long long seed)
{
    //Block index
    int bx = blockIdx.x;

    //Thread index
    int tx = threadIdx.x;

    //First element processed by the block
    int begin = NUM_THREADS_PER_BLOCK * bx;

	float s;	//step length

	unsigned int w = 0;
	float w_temp;

	PhotonStruct p = DeviceMem.p[begin+tx];
	Fibers f = DeviceMem.f[begin+tx];

	int new_layer;

	curandState state = DeviceMem.state[begin+tx];
    curand_init(seed, begin+tx, 0, &state);

	//First, make sure the thread (photon) is active
	unsigned int ii = 0;
	if(!DeviceMem.thread_active[begin+tx]) ii = NUMSTEPS_GPU;

	bool k = true;

	for(;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
	{
		if(layers_dc[p.layer].mutr!=FLT_MAX)
			s = -__logf(rnGen(&state))*layers_dc[p.layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
		else
			s = 100.0f;//temporary, say the step in glass is 100 cm.

		//Check for layer transitions and in case, calculate s
		new_layer = p.layer;
		if(p.z+s*p.dz<layers_dc[p.layer].z_min){new_layer--; s = __fdividef(layers_dc[p.layer].z_min-p.z,p.dz);} //Check for upwards reflection/transmission & calculate new s
		if(p.z+s*p.dz>layers_dc[p.layer].z_max){new_layer++; s = __fdividef(layers_dc[p.layer].z_max-p.z,p.dz);} //Check for downward reflection/transmission

		p.x += p.dx*s;
		p.y += p.dy*s;
		p.z += p.dz*s;

		if(p.z>layers_dc[p.layer].z_max)p.z=layers_dc[p.layer].z_max;//needed?
		if(p.z<layers_dc[p.layer].z_min)p.z=layers_dc[p.layer].z_min;//needed?

		if(new_layer!=p.layer)
		{
			// set the remaining step length to 0
			s = 0.0f;

			if(reflect(&p,new_layer,&state)==0u)//Check for reflection
			{
				if(new_layer == 0)
				{ //Diffuse reflectance
					detect(&p,&f);
					p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
				if(new_layer > *n_layers_dc)
				{	//Transmitted
					p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
			}
		}

		if(s > 0.0f)
		{
			if (layers_dc[p.layer].mua == 0.0f)
				p.Ab_z ++;
			// Drop weight (apparently only when the photon is scattered)
			//w_temp = __float2uint_rn(layers_dc[p.layer].mua*layers_dc[p.layer].mutr*__uint2float_rn(p.weight));
			w_temp = layers_dc[p.layer].mua*layers_dc[p.layer].mutr*p.weight;
			p.weight -= w_temp;
			spin(&p,DeviceG,&state);
		}

		if(p.Ab_z>10000)
			p.weight = 0;

		if(!photonSurvive(&p,&state)) //if the photon doesn't survive
		{
			k = false;
			if(atomicAdd(DeviceMem.num_terminated_photons,1u) < (*num_photons_dc-NUM_THREADS))
			{	// Ok to launch another photon
				launchPhoton(&p,&state);//Launch a new photon
			}
			else
			{	// No more photons should be launched.
				DeviceMem.thread_active[begin+tx] = 0u; // Set thread to inactive
				ii = NUMSTEPS_GPU;				// Exit main loop
			}
		}

	}//end main for loop!

	__syncthreads();//necessary?

	//save the state of the MC simulation in global memory before exiting
	DeviceMem.p[begin+tx] = p;	//This one is incoherent!!!
	DeviceMem.f[begin+tx] = f;

}//end MCd

__device__ void launchPhoton(PhotonStruct* p, curandState *state)
{
	float rnd_position, rnd_Azimuth, rnd_direction, rnd_rotated;
	float AzimuthAngle;
	float launchPosition;
	float theta_direction;
	float rotated_angle;
	float uxprime, uyprime, uzprime;
	float angle = -ANGLE * PI / 180;


	rnd_position   = rnGen(state);
	rnd_Azimuth    = rnGen(state);
	rnd_direction  = rnGen(state);
	rnd_rotated    = rnGen(state);
	AzimuthAngle   = 2 * PI * rnd_Azimuth;
	rotated_angle  = 2 * PI * rnd_rotated;

	float beam_width = 0.01;  // 100 um, Gaussian beam profile

	launchPosition = beam_width*sqrt(-log(rnGen(state))/2.0);

	p->x = launchPosition*cos(AzimuthAngle)/cos(angle);
	p->y = launchPosition*sin(AzimuthAngle);
	p->z = 0.0;

	theta_direction = asin(SOURCE_NA/N_SOURCE)*rnd_direction;
	p->dz = cos(theta_direction);
	p->dx = sin(theta_direction) * cos(rotated_angle);
	p->dy = sin(theta_direction) * sin(rotated_angle);

	uxprime = cos(angle)*p->dx - sin(angle)*p->dz;
	uyprime = sin(theta_direction)*sin(rotated_angle);
	uzprime = sin(angle)*p->dx + cos(angle)*p->dz;

	p->dx = uxprime, p->dy = uyprime, p->dz = uzprime;

	p->layer = 1;
	p->weight = *start_weight_dc; //specular reflection!
	p->Ab_z = 0;
}


__global__ void launchPhoton_Global(MemStruct DeviceMem, unsigned long long seed)
{
	int bx=blockIdx.x;
    int tx=threadIdx.x;

    //First element processed by the block
    int begin=NUM_THREADS_PER_BLOCK*bx;

	PhotonStruct p;

	curandState state = DeviceMem.state[begin+tx];
    curand_init(seed, 0, 0, &state);

	launchPhoton(&p,&state);

	//__syncthreads();//necessary?
	DeviceMem.p[begin+tx]=p;//incoherent!?
}

__device__ int binarySearch(float *data, float value)
{
    int middle;
	int left = 0, right = 180;
    while (left <= right)
    {
        middle = (right + left) / 2;

        if (data[middle] == value)
            return middle;

        if (data[middle] > value)
            right = middle - 1;
        else
            left = middle + 1;
    }
	if (data[middle] > value)
	    return middle;
	else
		return middle + 1;
}

__device__ void spin(PhotonStruct* p, G_Array g, curandState *state)
{
	float theta, cost, sint;	// cosine and sine of the
						// polar deflection angle theta.
	float cosp, sinp;	// cosine and sine of the
						// azimuthal angle psi.
	float temp;
	float tempdir=p->dx;

	float rn = rnGen(state);
    int sample;

	sample = binarySearch(g.cumf,rn);

	theta = sample-1+__fdividef((rn-g.cumf[sample-1]),(g.cumf[sample]-g.cumf[sample-1]));
  theta = __fdividef(theta*PI,180);
	cost = cos(theta);

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rnGen(state),&sinp,&cosp);// spin psi [0-2*PI)

	temp = sqrtf(1.0f - p->dz*p->dz);

	if(temp==0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost,p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp),temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp),temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(p->dx*p->dx+p->dy*p->dy+p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin



__device__ unsigned int reflect(PhotonStruct* p, int new_layer, curandState *state)
{
	//Calculates whether the photon is reflected (returns 1) or not (returns 0)
	// Reflect() will also update the current photon layer (after transmission) and photon direction (both transmission and reflection)

	float n1 = layers_dc[p->layer].n;
	float n2 = layers_dc[new_layer].n;
	float r;
	float cos_angle_i = fabsf(p->dz);

	if(n1==n2)//refraction index matching automatic transmission and no direction change
	{
		p->layer = new_layer;
		return 0u;
	}

	if(n1>n2 && n2*n2<n1*n1*(1-cos_angle_i*cos_angle_i))//total internal reflection, no layer change but z-direction mirroring
	{
		p->dz *= -1.0f;
		return 1u;
	}

	if(cos_angle_i==1.0f)//normal incident
	{
		r = __fdividef((n1-n2),(n1+n2));
		if(rnGen(state)<=r*r)
		{
			//reflection, no layer change but z-direction mirroring
			p->dz *= -1.0f;
			return 1u;
		}
		else
		{	//transmission, no direction change but layer change
			p->layer = new_layer;
			return 0u;
		}
	}

	//gives almost exactly the same results as the old MCML way of doing the calculation but does it slightly faster
	// save a few multiplications, calculate cos_angle_i^2;
	float e = __fdividef(n1*n1,n2*n2)*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
	r=2*sqrtf((1.0f-cos_angle_i*cos_angle_i)*(1.0f-e)*e*cos_angle_i*cos_angle_i);//use r as a temporary variable
	e=e+(cos_angle_i*cos_angle_i)*(1.0f-2.0f*e);//Update the value of e
	r = e*__fdividef((1.0f-e-r),((1.0f-e+r)*(e+r)));//Calculate r

	if(rnGen(state)<=r)
	{
		// Reflection, mirror z-direction!
		p->dz *= -1.0f;
		return 1u;
	}
	else
	{
		// Transmission, update layer and direction
		r = __fdividef(n1,n2);
		e = r*r*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
		p->dx *= r;
		p->dy *= r;
		p->dz = copysignf(sqrtf(1-e) ,p->dz);
		p->layer = new_layer;
		return 0u;
	}

}

__device__ unsigned int photonSurvive(PhotonStruct* p, curandState *state)
{
	//Calculate wether the photon survives (returns 1) or dies (returns 0)

	if(p->weight>WEIGHTI) return 1u; // No roulette needed
	if(p->weight==0.0f) return 0u;	// Photon has exited slab, i.e. kill the photon

	if(rnGen(state) < CHANCE)
	{
		//p->weight = __float2uint_rn(__fdividef((float)p->weight,CHANCE));
		p->weight = __fdividef(p->weight,CHANCE);
		return 1u;
	}
	return 0u;
}

__device__ void detect(PhotonStruct* p, Fibers* f)
{
	float angle = ANGLE*PI/180;
	float critical = asin(f->NA[1]/ N_DETECTOR);
    float uz_rotated=(p->dx*sin(angle))+(p->dz*cos(angle));
	float uz_angle = acos(fabs(uz_rotated));
	float distance;

	if(uz_angle <= critical)  // successfully detected
	{
		if(NORMAL)
		{
			//
			/*
			for(int i = 1; i <= 3 ; i++)
			{
				if(pow((p->x-f->position[i])*cos(angle),2) + pow(p->y,2) <= f->radius[i]*f->radius[i])
					f->data[i] += p->weight;
			}

			for(int i = 4; i <= 6 ; i++)
			{
				if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
					f->data[i] += p->weight;
			}
			*/

			// ISS annular
			distance = sqrt(p->x * p->x + p->y * p->y);

			for(int i = 1; i <= 6 ; i++)
			{
				if((distance>=(f->position[i]-f->radius[i])) && (distance<=(f->position[i]+f->radius[i])))
				{
					float temp;
					temp = (distance*distance + f->position[i]*f->position[i] - f->radius[i]*f->radius[i])/(2*distance*f->position[i]);
					// check for rounding error!
					if(temp > 1.0f)
						temp = 1.0f;

					f->data[i] += p->weight * acos(temp) * RPI;

				}
			}

		}
		else
		{
			for(int i = 1; i <= 3 ; i++)
			{
				if(pow((p->x-f->position[i])*cos(angle),2) + pow(p->y,2) <= f->radius[i]*f->radius[i])
					f->data[i] += p->weight;
			}

			for(int i = 4; i <= 6 ; i++)
			{
				if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
					f->data[i] += p->weight/2;
			}

			for(int i = 7; i <= 9 ; i++)
			{
				if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
					f->data[i] += p->weight/2;
			}
		}
	}
     return;
}

int initDCMem(SimulationStruct* sim)
{
	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(n_layers_dc,&(sim->n_layers),sizeof(unsigned int));

	// Copy start_weight_dc to constant device memory
	cudaMemcpyToSymbol(start_weight_dc,&(sim->start_weight),sizeof(float));

	// Copy layer data to constant device memory
	cudaMemcpyToSymbol(layers_dc,sim->layers,(sim->n_layers+2)*sizeof(LayerStruct));

	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(num_photons_dc,&(sim->number_of_photons),sizeof(unsigned int));

	return 0;
}

int initG(G_Array* HostG, G_Array* DeviceG, int index)
{
	HostG->all = (float*) malloc(181*sizeof(float));
	HostG->cumf = (float*) malloc(181*sizeof(float));
	HostG->prob = (float*) malloc(181*sizeof(float));
	cudaMalloc((void**)&DeviceG->all,181*sizeof(float));
	cudaMalloc((void**)&DeviceG->cumf,181*sizeof(float));
	cudaMalloc((void**)&DeviceG->prob,181*sizeof(float));
	inputG(index,HostG);
	cudaMemcpy(DeviceG->all,HostG->all,181*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(DeviceG->cumf,HostG->cumf,181*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(DeviceG->prob,HostG->prob,181*sizeof(float),cudaMemcpyHostToDevice);
	return 0;
}

int initMemStruct(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim)
{
	// Allocate p on the device!!
	cudaMalloc((void**)&DeviceMem->p,NUM_THREADS*sizeof(PhotonStruct));

	// Allocate thread_active on the device and host
	HostMem->thread_active = (unsigned int*) malloc(NUM_THREADS*sizeof(unsigned int));
	if(HostMem->thread_active==NULL){printf("Error allocating HostMem->thread_active"); exit (1);}
	for(int i=0;i<NUM_THREADS;i++)HostMem->thread_active[i]=1u;

	cudaMalloc((void**)&DeviceMem->thread_active,NUM_THREADS*sizeof(unsigned int));
	cudaMemcpy(DeviceMem->thread_active,HostMem->thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyHostToDevice);

	//Allocate num_launched_photons on the device and host
	HostMem->num_terminated_photons = (unsigned int*) malloc(sizeof(unsigned int));
	if(HostMem->num_terminated_photons==NULL){printf("Error allocating HostMem->num_terminated_photons"); exit (1);}
	*HostMem->num_terminated_photons=0;

	cudaMalloc((void**)&DeviceMem->num_terminated_photons,sizeof(unsigned int));
	cudaMemcpy(DeviceMem->num_terminated_photons,HostMem->num_terminated_photons,sizeof(unsigned int),cudaMemcpyHostToDevice);

	//Allocate and initialize fiber f on the device and host
	HostMem->f = (Fibers*) malloc(NUM_THREADS*sizeof(Fibers));
	cudaMalloc((void**)&DeviceMem->f,NUM_THREADS*sizeof(Fibers));
	initFiber(HostMem->f);
	cudaMemcpy(DeviceMem->f,HostMem->f,NUM_THREADS*sizeof(Fibers),cudaMemcpyHostToDevice);

	//Allocate states on the device and host
	cudaMalloc((void**)&DeviceMem->state,NUM_THREADS*sizeof(curandState));


	return 1;
}

void freeMemStruct(MemStruct* HostMem, MemStruct* DeviceMem)
{
	free(HostMem->thread_active);
	free(HostMem->num_terminated_photons);
	free(HostMem->f);

	cudaFree(DeviceMem->thread_active);
	cudaFree(DeviceMem->num_terminated_photons);
	cudaFree(DeviceMem->f);
	cudaFree(DeviceMem->state);
}

void freeG(G_Array* HostG, G_Array* DeviceG)
{
	free(HostG->all);
	free(HostG->cumf);
	free(HostG->prob);
	cudaFree(DeviceG->all);
	cudaFree(DeviceG->cumf);
	cudaFree(DeviceG->prob);
}
