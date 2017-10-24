#include "header.h"
#include <unistd.h>

static ofstream myfile;

int main(int argc,char* argv[])
{



	SimulationStruct* simulations;
	int n_simulations;
	unsigned long long seed = (unsigned long long) time(NULL);// Default, use time(NULL) as seed
	n_simulations = readData(&simulations); // read the input file

	if(n_simulations == 0)
	{
		printf("Something wrong with read_simulation_data!\n");
		return 1;
	}
	else
	{
		//printf("Successfully read data!\n");
	}

	clock_t time1,time2;

	// Start the clock
    time1=clock();




		char output_file[50] = "../output_DRS/simulation_0.txt";
		//CHECK IF file EXIST, IF YES, CREATE simulation_X+1.txt
		int file_index=1;
		while(access(output_file,F_OK)!= -1){
			memset(output_file,0,50);
			sprintf(output_file,"../output_DRS/simulation_%d.txt",file_index);
			file_index++;
		}
		myfile.open (output_file,ios::app);








	//perform all the simulations
	for(int i = 0; i < n_simulations; i++)
	{
		// Run a simulation
		printf("simulation #%d\n",i+1);
		doOneSimulation(&simulations[i],i,myfile);
	}

	myfile.close();
	printf("result saved in %s",output_file);

	time2=clock();
	printf("Simulation time: %.2f sec\n",(double)(time2-time1)/CLOCKS_PER_SEC);

	freeSimulationStruct(simulations, n_simulations);

	//system("PAUSE");
	return 0;
}

void freeSimulationStruct(SimulationStruct* sim, int n_simulations)
{
	for(int i = 0;i < n_simulations; i++) free(sim[i].layers);
	free(sim);
}
