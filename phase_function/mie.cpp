/*
 * We translated the Fortran code written by Tony Durkin at UT, Austin
 * into ANSI C. Sept. 14, 1993. See main() for authors' information.
 *
 * To compile, type: cc -o miesphr miesphr.c or acc -o miesphr miesphr.c
 * whichever compiler complies to ANSI C.
 */
/**********************************************************************
Program written to calculate the scattering coefficient (mus) and
scattering anisotropy (g) for wavelength region specified by the user
the program uses subroutine callbh written by Bohren et al:
	absorption and scattering of light by small particles.
**********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mie.h"



complex Cform(float rl, float im)
{
  complex     c;

  c.r = rl;
  c.i = im;
  return (c);
}

float Creal(complex C)
{
  return (C.r);
}

float Cimag(complex C)
{
  return (C.i);
}

float Cabs(complex C)
{
  return (sqrt(C.r * C.r + C.i * C.i));
}

complex Cadd(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r + C2.r;
  c.i = C1.i + C2.i;
  return (c);
}

complex Csub(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r - C2.r;
  c.i = C1.i - C2.i;
  return (c);
}

/* (a + ib)(c + id) = ac-bd + i(bc+ad) */
complex Cmulti(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r * C2.r - C1.i * C2.i;
  c.i = C1.r * C2.i + C1.i * C2.r;
  return (c);
}

/* (a + ib)/(c + id) = (a+ib)(c-id)/(c^2+d^2) */
complex Cdiv(complex C1, complex C2)
{
  float       temp;
  complex     c;

  temp = 1 / (C2.r * C2.r + C2.i * C2.i);
  c.r = (C1.r * C2.r + C1.i * C2.i) * temp;
  c.i = (C1.i * C2.r - C1.r * C2.i) * temp;
  return (c);
}

complex Cconj(complex C)
{
  complex     ctemp;

  ctemp.r = C.r;
  ctemp.i = -C.i;
  return (ctemp);
}

/**************************************************************************
subroutine BHMie calculates amplitude scattering matrix
elements and efficiencies for extinction, total scattering
and backscattering for a given size parameter and
relative refractive index
*************************************************************************/
void
BHMie(float *X,
      complex * RefRel,
      int *Nang,
      complex * S1,
      complex * S2,
      // Q -> efficiency
      float *Qext,
      float *Qsca,
      float *Qback,
      float *Ganiso)
{
  static complex cd[3000];
  static complex cy;
  static complex can, cbn;
  static complex cxi;
  static complex cxi0, cxi1;
  static complex can1, cbn1, can2, cbn2;
  complex     ctemp1, ctemp2, ctemp3;

  static float dang, apsi, ymod;
  static double qsca1;
  static float apsi0, apsi1;
  static int  j, n;
  static float p, t;
  static float theta[100];
  static int  nstop;
  static float xstop;
  static double dn;
  static float fn;
  static int  jj;
  static float pi[100];
  static float pi0[100], pi1[100];
  static double dx;
  static int  nn;
  static float rn;
  static float ganisotmp;
  static int  rn1;
  static float chi, mu[100], tau[100];
  static double psi;
  static int  nmx;
  static float chi0, chi1;
  static double psi0, psi1;

  dx = *X;
  cy = Cform((*X) * RefRel->r, (*X) * RefRel->i);

  /*************************************************************************
  series terminated after nstop terms
  *************************************************************************/
  xstop = *X + 4 * pow(*X, 0.3333) + 2.0;
  nstop = xstop;
  ymod = Cabs(cy);
  nmx = (xstop > ymod) ? xstop : ymod + 14;
  dang = PI / (2 * (*Nang - 1));

  for (j = 0; j < *Nang; j++) {
    theta[j] = j * dang;
    mu[j] = cos(theta[j]);
  }

  /************************************************************************
  logarithmic derivative cd(j) calculated by downward
  recurrence beginning with initial value 0.0+i*0.0
  at j=nmx
  ***********************************************************************/
  cd[nmx] = Cform(0.0, 0.0);
  nn = nmx;

  for (n = 0; n < nn; n++) {
    rn = nmx - n + 1;
    ctemp1 = Cform(rn, 0.0);
    ctemp1 = Cdiv(ctemp1, cy);
    ctemp2 = Cadd(cd[nmx - n], ctemp1);
    ctemp3 = Cform(1, 0);
    ctemp3 = Cdiv(ctemp3, ctemp2);
    cd[nmx - n - 1] = Csub(ctemp1, ctemp3);
  }

  for (j = 0; j < *Nang; j++) {
    pi0[j] = 0.0;
    pi1[j] = 1.0;
  }

  nn = (*Nang << 1) - 1;
  for (j = 0; j < nn; j++)
    S1[j] = S2[j] = Cform(0.0, 0.0);

  /************************************************************************
  riccati-bessel functions with real argument X
  calculated by upward recurrence
  *********************************************************************/
  psi0 = cos(dx);
  psi1 = sin(dx);
  chi0 = -sin(*X);
  chi1 = cos(*X);
  apsi0 = psi0;
  apsi1 = psi1;
  cxi0 = Cform(apsi0, -chi0);
  cxi1 = Cform(apsi1, -chi1);
  *Qsca = 0.0;
  *Ganiso = 0.0;

  can1 = Cform(0.0, 0.0);
  cbn1 = Cform(0.0, 0.0);
  can2 = Cform(0.0, 0.0);
  cbn2 = Cform(0.0, 0.0);
  qsca1 = 0.0;

  n = 1;
  do {
    dn = (double) n;
    rn = (float) n;
    fn = (2. * rn + 1.) / (rn * (rn + 1.));
    psi = (2. * dn - 1.) * psi1 / dx - psi0;
    apsi = psi;
    chi = (2. * rn - 1.) * chi1 / (*X) - chi0;
    cxi = Cform(apsi, -chi);

    ctemp1 = Cdiv(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    can = Cdiv(ctemp2, ctemp3);

    ctemp1 = Cmulti(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    cbn = Cdiv(ctemp2, ctemp3);

    can1 = can;
    cbn1 = cbn;
    *Qsca += (2 * rn + 1) * (Cabs(can) * Cabs(can) + Cabs(cbn) * Cabs(cbn));

    for (j = 0; j < *Nang; j++) {
      jj = (*Nang << 1) - j - 2;
      pi[j] = pi1[j];
      tau[j] = rn * mu[j] * pi[j] - (rn + 1) * pi0[j];
      p = pow(-1, n - 1);

      ctemp1 = Cform(fn * (can.r * pi[j] + cbn.r * tau[j]),
		     fn * (can.i * pi[j] + cbn.i * tau[j]));
      S1[j] = Cadd(S1[j], ctemp1);
      t = pow(-1., n);

      ctemp1 = Cform(fn * (can.r * tau[j] + cbn.r * pi[j]),
		     fn * (can.i * tau[j] + cbn.i * pi[j]));
      S2[j] = Cadd(S2[j], ctemp1);

      if (j == jj)
	continue;

      ctemp1 = Cform(fn * (can.r * pi[j] * p + cbn.r * tau[j] * t),
		     fn * (can.i * pi[j] * p + cbn.i * tau[j] * t));
      S1[jj] = Cadd(S1[jj], ctemp1);

      ctemp1 = Cform(fn * (can.r * tau[j] * t + cbn.r * pi[j] * p),
		     fn * (can.i * tau[j] * t + cbn.i * pi[j] * p));
      S2[jj] = Cadd(S2[jj], ctemp1);
    }

    psi0 = psi1;
    psi1 = psi;
    apsi1 = psi1;
    chi0 = chi1;
    chi1 = chi;
    cxi1 = Cform(apsi1, -chi1);
    rn1 = rn;
    n = n + 1;
    rn = n;

    for (j = 0; j < *Nang; j++) {
      pi1[j] = ((2. * rn - 1.) / (rn - 1.)) * mu[j] * pi[j];
      pi1[j] = pi1[j] - rn * pi0[j] / (rn - 1.);
      pi0[j] = pi[j];
    }

    dn = n;
    rn = n;
    fn = (2. * rn + 1.) / (rn * (rn + 1.));
    psi = (2. * dn - 1.) * psi1 / dx - psi0;
    apsi = psi;
    chi = (2. * rn - 1.) * chi1 / (*X) - chi0;
    cxi = Cform(apsi, -chi);

    ctemp1 = Cdiv(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    can = Cdiv(ctemp2, ctemp3);

    ctemp1 = Cmulti(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    cbn = Cdiv(ctemp2, ctemp3);

    can2 = can;
    cbn2 = cbn;

    ctemp1 = Cmulti(can1, Cconj(can2));
    ctemp2 = Cmulti(cbn1, Cconj(cbn2));
    ctemp3 = Cadd(ctemp1, ctemp2);
    ganisotmp = rn1 * (rn1 + 2.0) * Creal(ctemp3);

    *Ganiso += ganisotmp / (rn1 + 1.);
    ctemp1 = Cmulti(can1, Cconj(cbn1));
    *Ganiso += (2. * rn1 + 1.) * Creal(ctemp1) / (rn1 * (rn1 + 1.0));
  } while (n - 1 - nstop < 0);

  *Qsca = (2. / (*X * *X)) * *Qsca;
  *Qext = (4. / (*X * *X)) * Creal(S1[0]);
  *Qback = (4. / (*X * *X)) *
    Cabs(S1[2 * (*Nang) - 2]) * Cabs(S1[2 * (*Nang) - 2]);
  *Ganiso = (4. / (*X * *X)) * *Ganiso;
  *Ganiso = *Ganiso / *Qsca;
}

/**********************************************************************
This function first calculates the size parameter(x) and relative
refractive index(refrel) for a given sphere refractive index, medium
refractive index, radius and free space wavelength. It then calls BHMie
to compute amplitude scattering matrix elements and efficiencies.
**********************************************************************/
static complex s1[181], s2[181];
void
CallBH(float *Wavelen,
       float *Qsca,
       float *RefMed,		/* n of the surrounding medium. */
       float *RefRe,		/* n of the sphere [real part]. */
       float *Radius,
       float *Ganiso,
       float *Qext)
{
  static complex refrel;	/* relative n of the sphere [complex]. */
  static complex ctemp;

  static int  nang, nan, j;
  static float ang, dang;
  static float qback;
  static float x, aj, pol, s11nor, s11, s12, s33, s34;
  float       temp1, temp2;

  refrel = Cform(*RefRe / (*RefMed), 0.0);
  x = 2. * PI * (*Radius) * (*RefMed) / (*Wavelen);

  /******************************************************************
  nang=number of angles between 0 and 90 degrees
  matrix elements calculated at 2*nang-1 angles
  including 0, 90, and 180 degrees
  ******************************************************************/
  nang = 91;
  dang = PI / (2 * (nang - 1));
  BHMie(&x, &refrel, &nang, s1, s2, Qext, Qsca, &qback, Ganiso);


#if 0
  /******************************************************************
  s33 and s34 matrix elements normalized by s11
  s11 is normalized to 1.0 in the forward direction
  pol=degree of polarization(incident unpolirized light)
  ******************************************************************/
  temp1 = Cabs(s1[0]);
  temp2 = Cabs(s2[0]);
  s11nor = 0.5 * (temp1 * temp1 + temp2 * temp2);
  nan = (nang << 1) - 1;

  for (j = 0; j < nan; j++) {
    aj = j + 1;
    temp1 = Cabs(s1[j]);
    temp2 = Cabs(s2[j]);
    s11 = 0.5 * (temp2 * temp2 + temp1 * temp1);
    s12 = 0.5 * (temp2 * temp2 - temp1 * temp1);

    pol = -s12 / s11;
    ctemp = Cmulti(s2[j], Cconj(s1[j]));
    s33 = Creal(ctemp);
    s33 = s33 / s11;
    s34 = Cimag(ctemp);
    s34 = s34 / s11;
    s11 = s11 / s11nor;
    ang = dang * (aj - 1.) * 57.2958;
    /* printf("%14.6g %14.6g %14.6g %14.6g %14.6g", ang,s11,pol,s33,s34); */
  }
#endif
}

void MIE(MieStruct mie,int NUM_WAVE)
{
  double A = 1.5725;
  double B = 0.0031080;
  double C = 0.00034779;
  float       qext;		/* extinction efficiency. */
  float       number_density;	/* number density. [1/um^3] */

  float       qsca[NUM_WAVE];	/* scattering efficiency. */
  float       g[NUM_WAVE];	/* anisotropy, g. */
  int         i;


  //input through MieStruct
  float       rad;			/* radius of the sphere. [um] */
	float       specific_weight_spheres;/* [g/cc]. */
	float       specific_weight_solvent;/* [g/cc]. */
	float       concentration_by_weight;/* [g/g]. */
  float       wlen[NUM_WAVE];	/* wavelength. [um] */
  float       refre[NUM_WAVE];		/* refractive index of the sphere. */
  float       refre_water[NUM_WAVE];
  float       Hb[NUM_WAVE];

  FILE        *filep;


  //initialize
	rad = mie.rad;
	specific_weight_spheres = mie.specific_weight_spheres;
	specific_weight_solvent = mie.specific_weight_solvent;
	//concentration_by_weight = mie.concentration_by_weight;


  //initialize phantom's refractance and wavelength(wlen)
  for(int i = 0;i<NUM_WAVE;i++)
  {
    wlen[i] = mie.wave_length[i];
    //use refre = n(λ) = A + B/λ^2 + C/λ^2
    refre[i] = A + B/(pow(wlen[i]/1000,2)) + C/(pow(wlen[i]/1000,4));
  }
  //initialize water's refractance
  for(int i=0;i<NUM_WAVE;i++)
  {
    double temp = pow((wlen[i]*0.001),2);
    refre_water[i] = sqrt(0.566695982*temp/(temp-0.0055084151894)+0.1731900098*temp/(temp-0.01818488474)+0.02095951857*temp/(temp-0.02625439472)+0.1125228406*temp/(temp-10.73842352)+1);
  }
  //initialize Hb concentration
  char* fname = "mua_hb.txt";
  if (!(filep = fopen(fname, "r"))) {
    puts("Can't open file mua_hb.");
    exit(1);
  }
  int num_hb=800;
  float wl[num_hb],hb[num_hb];
  for(int i=0;i<num_hb;i++) fscanf(filep, "%f\t%f",&(wl[i]),&(hb[i]));
  for(int i=0;i<NUM_WAVE;i++)
  {
      for(int j=0;j<num_hb;j++)
      {
          if(wl[j]==wlen[i])
          {
              Hb[i] = hb[j];
              //for checking
              //printf("wave:%.0f, Hb:%f\n",wlen[i],Hb[i]);
              break;
          }
      }
  }


  fclose(filep);
















  rad /= 2;
  number_density = mie.num_density;
  //concentration_by_weight * specific_weight_solvent/ (specific_weight_spheres * 4 * PI / 3 * pow(rad, 3));


  for (i = 0; i < NUM_WAVE; i++) {
    //printf("%f\n",mie.wave_length[i]);
    wlen[i] /= 1000;
    CallBH(&(wlen[i]), &(qsca[i]), &(refre_water[i]), &(refre[i]), &rad, &(g[i]), &qext);

    //========output phase function====================

    float p[181];
    for(int j =0;j<181;j++){

      float temp1 = Cabs(s1[j]);
      float temp2 = Cabs(s2[j]);
      p[j] = 0.5 * (temp2 * temp2 + temp1 * temp1);

    }
    char fname_phase[32] = {};
    sprintf(fname_phase,"phase/phaseOutput_%.0fnm.txt",mie.wave_length[i]);
    if (!(filep = fopen(fname_phase, "w"))) {
      puts("Can't open phase/.");
      exit(1);
    }

    for (int j = 0; j < 181; j++)
	   fprintf(filep, "%14.6f\n",p[j]/100.0);

    fclose(filep);



  }
  //=========output origin mie.cpp======================================
  /* output results. */

    char* fname_standard = "input_data/mieOutput_standard.txt";
    if (!(filep = fopen(fname_standard, "w"))) {
      puts("Can't open file input_data/.");
      exit(1);
    }
    fprintf(filep, "%14s %14s %14s %14s\n",
	      "wavelen [um]", "Qsca", "mus [1/cm]", "g");

    for (i = 0; i < NUM_WAVE; i++)
	  fprintf(filep, "%14.6f %14.6f %14.6f %14.6f\n",
		wlen[i], qsca[i],
		qsca[i] * PI * rad * rad * number_density * 1E4, g[i]);

    fclose(filep);


//=======output index.txt=============================================

    char* fname3 = "input_data/index.txt";
    if (!(filep = fopen(fname3, "w"))) {
      puts("Can't open input_data/.");
      exit(1);
    }


    for (i = 0; i < NUM_WAVE; i++)
      fprintf(filep, "%14d %14.6f %14.6f\n",
      int(1000*wlen[i]),refre_water[i], refre[i]);

    fclose(filep);

//======output mua mus================================================

    char* fname4 = "input_data/input.txt";
    if (!(filep = fopen(fname4, "w"))) {
      puts("Can't open file.");
      exit(1);
    }

    fprintf(filep, "1.00\n");//set up thickness
    for (i = 0; i < NUM_WAVE; i++){
      float mus = qsca[i] * PI * rad * rad * number_density * 1E4; //1/cm
      float mua = Hb[i];
      fprintf(filep, "%.6f %14.6f %14.6f %14.6f\n",mua,mus,mua,mus);
    }
    fclose(filep);


}
