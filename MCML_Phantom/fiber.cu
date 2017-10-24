#include "header.h"

void initFiber(Fibers* f)
{
	for(int i = 0; i < NUM_THREADS; i++)
	{
	    for(int j = 0; j <= NUM_DETECTOR; j++)
			f[i].data[j] = 0;

		f[i].radius[0]   = ILLUMINATION_R;          // source fiber
	    f[i].NA[0]       = SOURCE_NA;
	    f[i].angle[0]    = ANGLE*PI/180;
	    f[i].position[0] = 0.0;

		if(NORMAL)
		{	/*
			f[i].radius[1]   = COLLECT_R;
			f[i].NA[1]       = DETECTOR_NA;
			f[i].position[1] = 0.024;
			f[i].angle[1]    = ANGLE*PI/180;

			f[i].radius[2]   = COLLECT_R;
			f[i].NA[2]       = DETECTOR_NA;
			f[i].position[2] = 0.049;
			f[i].angle[2]    = ANGLE*PI/180;

			f[i].radius[3]   = COLLECT_R;
			f[i].NA[3]       = DETECTOR_NA;
			f[i].position[3] = 0.076;
			f[i].angle[3]    = ANGLE*PI/180;

			f[i].radius[4]   = COLLECT_R;
			f[i].NA[4]       = DETECTOR_NA;
			f[i].position[4] = 0.026;
			f[i].angle[1]    = ANGLE*PI/180;

			f[i].radius[5]   = COLLECT_R;
			f[i].NA[5]       = DETECTOR_NA;
			f[i].position[5] = 0.054;
			f[i].angle[5]    = ANGLE*PI/180;

			f[i].radius[6]   = COLLECT_R;
			f[i].NA[6]       = DETECTOR_NA;
			f[i].position[6] = 0.078;
			f[i].angle[6]    = ANGLE*PI/180;
			*/

			f[i].radius[1]   = COLLECT_R;
			f[i].NA[1]       = DETECTOR_NA;
			f[i].position[1] = 0.022;//0.25 (eso.) unit [cm]
			f[i].angle[1]    = ANGLE*PI/180;

			f[i].radius[2]   = COLLECT_R;
			f[i].NA[2]       = DETECTOR_NA;
			f[i].position[2] = 0.042;//0.51 (eso.) unit [cm]
			f[i].angle[2]    = ANGLE*PI/180;

			f[i].radius[3]   = COLLECT_R;
			f[i].NA[3]       = DETECTOR_NA;
			f[i].position[3] = 0.065;//0.82 (eso.) unit [cm]
			f[i].angle[3]    = ANGLE*PI/180;

			f[i].radius[4]   = COLLECT_R;
			f[i].NA[4]       = DETECTOR_NA;
			f[i].position[4] = 0.022;//0.25 (eso.) unit [cm]
			f[i].angle[1]    = ANGLE*PI/180;

			f[i].radius[5]   = COLLECT_R;
			f[i].NA[5]       = DETECTOR_NA;
			f[i].position[5] = 0.042;//0.51 (eso.) unit [cm]
			f[i].angle[5]    = ANGLE*PI/180;

			f[i].radius[6]   = COLLECT_R;
			f[i].NA[6]       = DETECTOR_NA;
			f[i].position[6] = 0.065;//0.82 (eso.) unit [cm]
			f[i].angle[6]    = ANGLE*PI/180;

		}
		else
		{
			for(int j = 1; j <= 3; j++)
			{
				f[i].radius[j]   = COLLECT_R;
				f[i].NA[j]       = DETECTOR_NA;
				f[i].position[j] = 0.032*j;
				f[i].angle[j]    = ANGLE*PI/180;
			}


			for(int j = 4; j <= 6; j++)
			{
				f[i].radius[j]   = COLLECT_R;
				f[i].NA[j]       = DETECTOR_NA;
				f[i].position[j] = 0.022*(j-3);
				f[i].angle[j]    = ANGLE*PI/180;
			}

			for(int j = 7; j <= 9; j++)
			{
				f[i].radius[j]   = COLLECT_R;
				f[i].NA[j]       = DETECTOR_NA;
				f[i].position[j] = -0.022*(j-6);
				f[i].angle[j]    = ANGLE*PI/180;
			}
		}

	}
}
