
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

struct ALFAOFF {
  double az;
  double za;
};

static struct ALFAOFF *alfaoff = NULL;

double roty( double, double, double);
double rotx( double, double, double);
int find_alfarot(double, int, double *, double * );
double find_alfaoff(int, int);
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

double find_alfaoff(beam, isza)
int beam;
int isza;
{
  int count;
  char *s;
  char *delim = " \t\n";
  FILE *fd;
  char line[80];

  if( !alfaoff ) {
    if( !(fd = fopen("/home/cole/deneva/alfasplit_AO/alfa_offsets.tcl", "r" )))
      {
	printf("FIND_ALFAOFF: no alfa_offsets.tcl file, returning\n");
	return(0.0);
      }

    alfaoff = (struct ALFAOFF *)malloc( sizeof(struct ALFAOFF)*7);
    bzero( alfaoff, sizeof(struct ALFAOFF)*7);
    count = 0;
    while( fgets( line, 80, fd )) {
      s = strtok( line, delim );
      if( !s || strcmp( s, "set" ) != 0 )
        continue;
      s = strtok( NULL, delim );
      if( !s || strcmp( s, "point(alfa_offs_table)" ) != 0 )
        continue;
      while( fgets( line, 80, fd )) {
        if(!(s = strtok( line, delim)))
          continue;
        if(!(s = strtok( NULL, delim)))
          continue;
        alfaoff[count].az = atof(s)/3600.0;
        if(!(s = strtok( NULL, delim)))
          continue;
        alfaoff[count].za = atof(s)/3600.0;
        count++;
        if( count >= 7 )
          break;
      }
      if( count >= 7 )
        break;
    }
  }


  if ( beam < 0 || beam > 6 )
    return 0.0;
  else if( isza )
    return alfaoff[beam].za;
  else
    return alfaoff[beam].az;
}

find_alfarot(angle, beam, paz, pza )
double angle;
int beam;
double *paz;
double *pza;
{
  double az, za;

  az  = find_alfaoff(beam, 0);
  za  = find_alfaoff(beam, 1);

 *paz = rotx( az, za, angle );
 *pza = roty( az, za, angle );
}

#define ALFARATIO (329.06/384.005)

/* negative angle associated with clockwize rotation */

double rotx( offx, offy, ang )
double offx, offy, ang;
{
  double ret, h, a1;

  offy *= ALFARATIO;

  h = sqrt( offx*offx + offy*offy );
  a1 = atan2( offy, offx);
  ret = h * cos( a1 - ang/180.0*M_PI );


  return ret;
}

/* adjust for the ellipse */
/* negative angle associated with clockwize rotation */

double roty( offx, offy, ang )
double offx, offy, ang;
{
  double ret, h, a1;

  offy *= ALFARATIO;
  h = sqrt( offx*offx + offy*offy );
  a1 = atan2( offy, offx);
  ret = h * sin( a1 - ang/180.0*M_PI );
  ret /= ALFARATIO;
  return ret;
}

/* given 
  J2000 ra/dec in hh.hhh ddd.ddd
  lst time in hours from midnight
  current epoch as in 2004.56789
  rotation angle
  alfa beam

  compute new ra/dec based on alfa beam offsets

*/

alfa_position( ra, dec, lst, epoch, angle, off1, off2, beam, pra, pdec )
double ra;
double dec;
double lst;
double epoch;
double angle;
double off1;
double off2;
int beam;
double *pra;
double *pdec;
{
  double rra, rdec, nra, ndec, bra, bdec, offaz, offza;
  double saz, sza, nsaz, nsza;

  rra = ra*M_PI*2.0/24.0;
  rdec = dec*M_PI*2.0/360.0;
  sla_preces("FK5",2000.0,epoch, &rra, &rdec );
  nra = rra * 24.0 / 2.0 / M_PI;
  ndec = rdec * 180.0 / M_PI;

  compute_azza( nra, ndec, lst, &saz, &sza );

  find_alfarot(angle, beam, &offaz, &offza );

  saz = saz - (offaz - off1)/deg_sin(sza);
  sza = sza - (offza - off2);
  compute_radec( saz, sza, lst, &bra, &bdec );

  rra = bra*2.0*M_PI/24.0;
  rdec = bdec*M_PI/180.0;
  sla_preces("FK5",epoch, 2000.0, &rra, &rdec );
  bra = rra * 24.0 / 2.0 / M_PI;
  bdec = rdec * 180.0 / M_PI;
  *pra = bra;
  *pdec = bdec;

}

