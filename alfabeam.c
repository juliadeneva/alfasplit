
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"

/*
  inputs

   ra  - hours hh.hhhh
   dec - degrees ddd.ddd
   lst - fractional hours hh.hhhhh
   mjd - fractional days
   rot - rotation angle in degrees.
   
*/
double find_alfaoff(int,int);
double rotx(double,double,double);
double roty(double,double,double);
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

main(argc, argv)
int argc;
char *argv[];
{
  double current_time, ra, dec, nra, ndec, saz, sza, lst, mjd, rra, rdec;
  double nsza, nsaz, bra, bdec, rot_ang, offx, offy;
  double offaz, offza;
  int beam, i;
  char *obs_date;
  int secs, test;

  ra = 0.0;
  dec = 0.0;
  lst = 0.0;
  mjd = 0.0;
  test = 0;
  rot_ang = 0;
  offaz = 0.0;
  offza = 0.0;
  for( i=1; i<argc; i++ ) {
    if ( strcmp( argv[i], "-ra" ) == 0 ) 
      ra = atof(argv[++i]);
    else if ( strcmp( argv[i], "-dec" ) == 0 ) 
      dec = atof(argv[++i]);
    else if ( strcmp( argv[i], "-lst" ) == 0 ) 
      lst = atof(argv[++i]);
    else if ( strcmp( argv[i], "-mjd" ) == 0 ) 
      mjd = atof(argv[++i]);
    else if ( strcmp( argv[i], "-rot" ) == 0 ) 
      rot_ang = atof(argv[++i]);
    else if ( strcmp( argv[i], "-offaz" ) == 0 ) 
      offaz = atof(argv[++i]);
    else if ( strcmp( argv[i], "-offza" ) == 0 ) 
      offza = atof(argv[++i]);
    else if ( strcmp( argv[i], "-test" ) == 0 ) 
      test = 1;
  }

  if( test ) {
    obs_date = "20040801";
    secs = 0;
    ra = 6.5;
    dec = 12.1;
    lst = 7.9;
    mjd = obs2mjd(obs_date);
    rot_ang = 0.0;
  } else {
    if( ra == 0.0 || dec == 0.0 || lst == 0.0 || mjd == 0.0 ) {
      fprintf( stderr, 
    "alfabeam -ra hh.hhh -dec ddd.ddd -lst hh.hhh -mjd ddd.ddddd -rot ang -offaz ddd.ddd -offza ddd.ddd\n" );
      exit(1);
    }
  }
    
/* inputs */

  current_time = 2000.0 + (mjd - obs2mjd("20000101") 
    + (secs/3600.0 - 4.0)/24.0)/365.25; 

  printf( "computing ALFA beam offsets for\n" );
  printf( "ra %f dec %f\n", ra, dec );
  printf( "lst %f mjd %f epoch_time %f\n", lst, mjd, current_time );
  printf("beam    ra(hr)            dec(deg)\n");

  for( beam=0; beam<7; beam++ ) {
    alfa_position( ra, dec, lst, current_time, rot_ang, offaz, offza, 
        beam, &bra, &bdec );
    printf("%d     %13.9f     %13.9f\n", beam, bra, bdec );
  }
}

/*  year:mm:dd */

int obs2mjd( dt )
char *dt;
{
  int year;
  int mon;
  int day;
  char date[80];

  strncpy( date, dt, 80 );

  date[8] = 0;
  day = atoi( &date[6] );
  date[6] = 0;
  mon = atoi( &date[4] );
  date[4] = 0;
  year = atoi( &date[0] );

  return( cal2mjd( year, mon, day ));
}

/*
  rra = ra*M_PI*2.0/24.0;
  rdec = dec*M_PI*2.0/360.0;
  sla_preces("FK5",2000.0,current_time, &rra, &rdec );
  nra = rra * 24.0 / 2.0 / M_PI;
  ndec = rdec * 180.0 / M_PI;


  compute_azza( nra, ndec, lst, &saz, &sza );

  for( beam=0; beam<7; beam++ ) {
    offx = find_alfaoff( beam, 0 );
    offy = find_alfaoff( beam, 1 );
    nsaz = saz + rotx( offx, offy, rot_ang );
    nsza = sza + roty( offx, offy, rot_ang );
    nsaz += offaz;
    nsza += offza;

    compute_radec( nsaz, nsza, lst, &bra, &bdec );

    rra = bra*2.0*M_PI/24.0;
    rdec = bdec*M_PI/180.0;
    sla_preces("FK5",current_time, 2000.0, &rra, &rdec );
    bra = rra * 24.0 / 2.0 / M_PI;
    bdec = rdec * 180.0 / M_PI;
*/
