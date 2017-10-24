#define PI                 3.141592654f


typedef struct
{
	float refmed;		/* refractive index of the surrouding medium. */
	float rad;			/* radius of the sphere. [um] */
	float specific_weight_spheres;/* [g/cc]. */
	float specific_weight_solvent;/* [g/cc]. */
	float concentration_by_weight;/* [g/g]. */
	float* wave_length;	/* an array store all the wavelength  */
	float num_density;
}MieStruct;

typedef struct {
  float       r, i;
}complex;

void MIE(MieStruct mie,int NUM_WAVE);

complex Cform(float rl, float im);
float Creal(complex C);
float Cimag(complex C);
float Cabs(complex C);


complex Cadd(complex C1, complex C2);
complex Csub(complex C1, complex C2);
complex Cmulti(complex C1, complex C2);
complex Cdiv(complex C1, complex C2);
complex Cconj(complex C);
void BHMie(float *X,
      complex * RefRel,
      int *Nang,
      complex * S1,
      complex * S2,
      // Q -> efficiency
      float *Qext,
      float *Qsca,
      float *Qback,
      float *Ganiso);
void CallBH(float *Wavelen,
             float *Qsca,
             float *RefMed,		/* n of the surrounding medium. */
             float *RefRe,		/* n of the sphere [real part]. */
             float *Radius,
             float *Ganiso,
             float *Qext);
