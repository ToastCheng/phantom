#include "mie.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char** argv){
        //read wavelength
        FILE *wave = fopen("wavelength.txt","r");
        if(!wave){
          perror("cannot open \"wavelength.txt\"\n");
        }
        int wavelength[100];
        int count=0;
        while(!feof(wave)){
          fscanf(wave,"%d\n",&(wavelength[count]));
          count++;
        }


        float size;
        float density;
        size = (argc>1) ? atof(argv[1]) : 1.0;
        density = (argc>2) ? atof(argv[2]) : 0.0204;
        printf("size:%f\ndensity:%f\n",size,density);
        MieStruct mie;
        //mie.refmed = 1.33; calculate respect to wavelength in mie.cpp
        mie.rad = size;
        mie.specific_weight_spheres = 1.05;
        mie.specific_weight_solvent = 1;
        //mie.concentration_by_weight = concentration;
        mie.num_density = density;
        //float wave_length[NUM_WAVE] = {400,410,414,420,430,440,450,460,470,480,
        //490,500,510,520,530,540,542,550,560,570,576,580,590,600,610,620,
        //630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800};
        //{700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000};
        mie.wave_length = (float*)malloc(sizeof(int)*count);
        for(int i=0;i<count;i++) mie.wave_length[i] = wavelength[i];
        MIE(mie,count);

        return 0;
}
