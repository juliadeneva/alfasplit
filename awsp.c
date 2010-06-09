/*--------------------------------------------------*
 * Important signal processing routines
 *
 *   Note - ALL of these functions work on raw correlation values
 *  which have been normalized such that  -1 <rho < 1
 *  i.e. any offset is removed and then it is divided by the
 *  maximum possible correlation.   
 *  This normalized value will (in general) be refered to as rho
 *  while the zero lag version will be zho
 *
 * Optimum thresholds (vthres = first threshold above 0)
 *    and resulting zerolag (normalized to 0:1)
 *   vthres 3 level = 0.6120003/sigma  - zho = 
 *   vthres 9 level = 0.266916/sigma -> zho = 0.2125855237
 *--------------------------------------------------*/
#include <stdio.h> 
#include <math.h> 

double inverse_cerf(double);
double attndb3lev(double);
double attndb9lev(double);
double pow3lev(double);
double pow9lev(double);
int vanvleck3lev(double *,int);
int vanvleck3levf(float *,int);
void vanvleck9lev(double *,int);
void vanvleck9levf(float *,int);
int makefactor2(int);
double hamming(int,int);


#ifndef M_PI
#define M_PI 3.1415926535
#endif

/*-----------------------------------------------------------------------*
 * Compute Necessary Adjustment to Power Levels for 3 level sampling
 *  
 * Note: Optimum Level Vthres = 0.612003 
 *     or 20*log10(0.6120003) = -4.26497 
 * This routine routines the required change in ATTENUATION required to achieve
 *   optimum performance.
 * THUS a positive result means MORE attenuation is required (reduce power)
 *      a negative result means LESS attenuation is required (increase power)
 *
 *-----------------------------------------------------------------------*/
double attndb3lev(double zho)
{
  double db;
  db = -4.26497 - 20*log10(sqrt(2)*inverse_cerf(zho));
  return(db);
}
/*-------------------------------------------------------------------*
 * Converts zerolag into an unbiased estimate of input power (sigma^2)
 * output is unitless: (sigma/vthre)^2
 * NOTE this value is scaled by vthres, such that pow =1 when threshold is
 * optimum
 *--------------------------------------------------------------------*/
double pow3lev(double zho)
{
  double tmp;
  tmp = inverse_cerf(zho);
  return(0.3745443672/(tmp*tmp*2));
}
/*-----------------------------------------------------------------------*
 * Convert zerolag into an unbiased estimate of power
 * NOTE this value is scaled by vthres, such that pow =1 when threshold is
 * optimum
 *-----------------------------------------------------------------------*/
double pow9lev(double zho)
{
  static double coef[8] = {-0.03241744594,4.939640303,-5.751574913, 34.83143031,-78.66637472, 213.7108496, -317.1011469, 245.8618017  };
  double pow;
  zho *= 0.0625;
  pow = ((((((coef[7]*zho + coef[6])*zho + coef[5])*zho +coef[4])*zho + coef[3])*zho + coef[2])*zho +coef[1])*zho +coef[0];
  return(pow);
}
/*-----------------------------------------------------------------------*
 * Compute Necessary Adjustment to Power Levels for 9 level sampling
 *   Optimum zerolag occurs when sigma = 1 and
 *      x1 = 0.266916 and vthres = +/- {1,3,5,7)*x1
 *      i.e. {0.266916, 0.800748, 1.33458, 1.86841} 
 *    see samp9lev.c for details
 *   this will produce a zho = 0.2125855237
 * 
 * Note a positive result means MORE attenuation is needed (power is high)
 *      a negative result means LESS attenuation is needed (power is low)
 *
 *-----------------------------------------------------------------------*/
double attndb9lev(double zho)
{
  static double coef[8] = {-0.03241744594,4.939640303,-5.751574913, 34.83143031,-78.66637472, 213.7108496, -317.1011469, 245.8618017  };
  double db,pow;
  zho *=0.0625;
  pow = ((((((coef[7]*zho + coef[6])*zho + coef[5])*zho +coef[4])*zho + coef[3])*zho + coef[2])*zho +coef[1])*zho +coef[0];
  db = 10*log10(pow);
  return(db);
}
/*----------------------------------------------------------------*
 *  round an integer UP to the next factor of two
 * this used to ensure an fft is applied on factors of 2, 
 * even if number of lags is an odd value
 *   Thus 1000 -> 1024 or 33 -> 64, etc
 *----------------------------------------------------------------*/
int makefactor2(int in)
{
   double tmp;
   int if2;
   tmp = ceil( log10((double)in)/log10(2.0));
   if2 = rint(pow(2.0,tmp));
   return(if2); 
}
/*-----------------------------------------------------------------*
 * Hamming Window Function - Standard time domain version:
 *              +  +
 *           +        + 
 *     +  +              +  +
 *   ----------------------------------------
 *     0  1  2  3  4  5  6  7   npts = 8
 *   Return the hanning window factor for index x, 
 * note index is "C" type i.e. index = 0,1,2, ... N-1
 *  Reference "Signals and Systems" by Ziemer, et al
 *-----------------------------------------------------------------*/
double hamming(int index, int npts)
{
  return( 0.54 - 0.46*cos(2*M_PI*((double)index + 0.5)/((double)npts)));
}
#define NO    0
#define YES   1

void *init_vanvleck(int);
int vanvleck_3lev(void *, double *, int );
int vanvleck_9lev(void *, double *, int );
int free_vanvleck(void *);

static void *vv_init = NULL;
static int vv_npts_save = 0;

int vanvleck3lev(double *rho,int npts)
{

  if( vv_npts_save != npts ) {
    if( vv_init )
      free_vanvleck( vv_init );
    vv_init = init_vanvleck(npts);
    vv_npts_save = npts;
  }
  vanvleck_3lev(init, rho, npts );
}

int vanvleck9lev(double *rho,int npts)
{
  if( npts_save != npts ) {
    if( init )
      free_vanvleck( init );
    init = init_vanvleck(npts);
    vv_npts_save = npts;
  }
  vanvleck_9lev(init, rho, npts );
}


int vanvleck3levf(float *rho,int npts)
{
  static double *dub = NULL;
  int i;

  if( npts_save != npts ) {
    if( init ) {
      free_vanvleck( init );
    }
    init = init_vanvleck(npts);

    if( dub )
      free(dub);

    dub = (double *)malloc(sizeof(double)*npts);
    vv_npts_save = npts;
  }
  vanvleck_3lev(init, dub, npts );
  for( i=0; i<npts; i++ ) 
    rho[i] = dub[i];
}

int vanvleck9levf(float *rho,int npts)
{
  static double *dub = NULL;
  int i;

  if( npts_save != npts ) {
    if( init ) {
      free_vanvleck( init );
    }
    init = init_vanvleck(npts);

    if( dub )
      free(dub);

    dub = (double *)malloc(sizeof(double)*npts);
    vv_npts_save = npts;
  }
  vanvleck_9lev(init, dub, npts );
  for( i=0; i<npts; i++ ) 
    rho[i] = dub[i];
}


/*---------------------------------------------------------------------*
 * return modified julian date given a gregorian calendar - based on slalib
 * routine sla_cldj which itself is based on Hatcher (1984) QJRAS 25, 53-55 
 * creation date 1999/07/10 [mjd=51639] (dunc@naic.edu) 
 *---------------------------------------------------------------------*/
double cal2mjd(int iy, int im, int id) 
{
  int leap;
  /* month lengths in days for normal and leap years */
  static int mtab[2][13] = {
    {0,31,28,31,30,31,30,31,31,30,31,30,31},
    {0,31,29,31,30,31,30,31,31,30,31,30,31}
  };

  /*validate year*/
  if (iy<-4699) {
    fprintf(stderr,"Invalid year passed to cal2mjd\n");
    return(0.0);
  } else {
    /* validate month */
    if (im<1 || im>12) {
      fprintf(stderr,"Invalid month passed to cal2mjd\n");
      return(0.0);
    } else {
      /* allow for leap year */
      leap = iy%4 == 0 && iy%100 != 00 || iy%400 == 0;
      /* validate day */
      if (id<1 || id>mtab[leap][im]) {
	fprintf(stderr,"Invalid day passed to cal2mjd\n");
        return(0.0);
      }
    }
  }
  return ((1461*(iy-(12-im)/10+4712))/4+(5+306*((im+9)%12))/10
           -(3*((iy-(12-im)/10+4900)/100))/4+id-2399904);
}

