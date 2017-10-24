#include "header.h"



void inputG(int index, G_Array *g)
{
	extern int wavelength[];

	//char *address = "C:\\Users\\user\\Desktop\\GPU\\Anglepattern\\1.0_20151001\\"; // can be modified

	char filename[100]; // complete address for the file
	ifstream infile;
	int i = wavelength[index];  // calculate which file to open
	sprintf(filename, "../phase_function/phase/phaseOutput_%dnm.txt", i);//sprintf(filename,"%s1.0_%dnm.txt",address,i);  // can be modified
    infile.open(filename);

	for(int i = 0; i < 181; i++)   // initialize all the g data
	{
	    g->all[i]=0.0;
		g->cumf[i]=0.0;
		g->prob[i]=0.0;
    }

    float sum = 0;

    for(int i = 0; i < 181; i++)
	{
		infile >> g->all[i];
        sum += g->all[i]*sin(i*PI/180);
    }

    for(int i = 0; i < 181; i++)
		g->prob[i] = g->all[i]*sin(i*PI/180)/sum;

    for(int i = 0; i < 181; i++)
	{
		if(i == 0)
		    g->cumf[i] = g->prob[i];
        else
			g->cumf[i] += g->cumf[i-1] + g->prob[i];
	}

    infile.close();
}

void outputFiber(SimulationStruct* sim, float *data, ofstream& myfile)
{


	float scale1 = (float)sim->number_of_photons;

	if(NORMAL)
	{
		for(int i = 0; i < 3; i++)
		{
			if (i==5) { myfile << data[i] / scale1; }
			else{ myfile << data[i] / scale1 << "\t"; }
		}
	}
	else
	{
		for(int i = 0; i < 3; i++)
		{
			if (i == 5) { myfile << data[i] / scale1; }
			else{ myfile << data[i] / scale1 << "\t"; }
		}
	}
	myfile << endl;

}


int readData(SimulationStruct** simulations)
{
	// parameters to be modified
	unsigned long number_of_photons = NUMBER_PHOTON ;
	const int n_simulations = NUMBER_SIMULATION;
	int n_layers = 2;                                   // double layer, default value = 2
	//float medium_n = 1.33;                            // refractive index of medium
	//float tissue_n = 1.60;                            // refractive index of tissue

	float start_weight;
	float upper_thickness;


	// read the file
	fstream myfile;
	myfile.open ("../phase_function/input_data/input.txt");
	float up_mua[n_simulations],up_mus[n_simulations],down_mua[n_simulations],down_mus[n_simulations];
	myfile >> upper_thickness;
	for(int i = 0; i < n_simulations; i++)
	    myfile >> up_mua[i] >> up_mus[i] >> down_mua[i] >> down_mus[i];
	myfile.close();

	fstream file;
	file.open("../phase_function/input_data/index.txt");
	float wavelength[n_simulations], tissue_n[n_simulations], medium_n[n_simulations];
	for(int i = 0; i < n_simulations; i++)
	    file >> wavelength[i] >> tissue_n[i] >> medium_n[i];

	file.close();

	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*) malloc(sizeof(SimulationStruct)*n_simulations);
	if(*simulations == NULL){perror("Failed to malloc simulations.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

	for(int i = 0;i < n_simulations; i++)
	{
		(*simulations)[i].number_of_photons=number_of_photons;
		(*simulations)[i].n_layers = n_layers;

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*) malloc(sizeof(LayerStruct)*(n_layers+2));
		if((*simulations)[i].layers == NULL){perror("Failed to malloc layers.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

		// Set upper refractive index (medium)
		(*simulations)[i].layers[0].n = medium_n[i];

		// Set the parameters of tissue (upper layer)
		(*simulations)[i].layers[1].n     = tissue_n[i];
		(*simulations)[i].layers[1].mua   = up_mua[i];
		(*simulations)[i].layers[1].z_min = 0;
		(*simulations)[i].layers[1].z_max = upper_thickness;
		(*simulations)[i].layers[1].mutr  = 1.0f/(up_mua[i]+up_mus[i]);

		// Set the parameters of tissue (lower layer)
		(*simulations)[i].layers[2].n     = tissue_n[i];
		(*simulations)[i].layers[2].mua   = down_mua[i];
		(*simulations)[i].layers[2].z_min = upper_thickness;
		(*simulations)[i].layers[2].z_max = 1.0;            // set as infinity
		(*simulations)[i].layers[2].mutr  = 1.0f/(down_mua[i]+down_mus[i]);

		// Set lower refractive index (medium)
		(*simulations)[i].layers[n_layers+1].n = tissue_n[i];

		//calculate start_weight
		//float n1=(*simulations)[i].layers[0].n;
		float n1 = N_SOURCE;
		float n2=(*simulations)[i].layers[1].n;
		float r = (n1-n2)/(n1+n2);
		r = r*r;
		//start_weight = (unsigned int)((double)0xffffffff*(1-r));
		start_weight = 1.0f-r;
		//printf("Start weight=%f\n",start_weight);
		(*simulations)[i].start_weight=start_weight;
	}
	return n_simulations;
}
