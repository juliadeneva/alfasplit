
/* create two files from one WAPP PULSAR file taken
   with ALFA.

   - use Desh defined filename
   - fix the header positions to the alfa beam

   Modified by JD, 20050331
   --get correct wapp # regardless of the filename
   --write ra & dec in split headers regardless of whether beam coordinates
   were written to the raw header
   --fixed beam indexing
*/

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "key.h"
int new_filename( char *, struct HEADERP *, char **);
void fix_positions( double, double, struct HEADERP *, int );
double deg2dms( double );
double dms2deg( double );
double find_alfaoff( int, int); 
double rotx(double,double,double);
double roty(double,double,double);
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();
int sla_preces(char *, double, double, double *, double *);
int compute_radec( double, double, double, double *, double *);
int compute_azza( double, double, double, double *, double *);
int alfa_position( double ra, double dec, double lst, double epoch, 
double angle, double off1, double off2, int beam, double *pra, double *pdec);

int new_filename(char *old, struct HEADERP *h, char *name[]);
int find_index(char *s);
int find_wappno(char *s);
int obs2mjd(char *dt);
void fix_positions(double ra, double dec, struct HEADERP *h, int beam);
double deg2dms(double t);
double dms2deg(double in);
int ast_seconds(char *s);
double rotx(double offx, double offy, double ang);
double roty(double offx, double offy, double ang);

int find_wapp(char* filename);

#define MAXFILENAME 8192

main(argc, argv )
int argc;
char *argv[];
{
  double ra, dec;
  char *name[2];
  struct HEADERP *h;
  unsigned char *data;
  int header_size, nifs, num_lags, datalen, i, wdatalen, rdatalen;
  int isalfa, off, lagformat, datasz, n, wapp;
  struct HEADERVAL header_val, isalfa_val;
  FILE *fd[2];
  double *beam_decj, *beam_raj, *beam_za, *beam_az;
  char frontend[16];

  if( !argv[1] ) {
    fprintf(stderr, "%s wappalfafile\n", argv[0] );
  }

  name[0] = ( char *)malloc( MAXFILENAME );
  name[1] = ( char *)malloc( MAXFILENAME );

  h = head_parse( argv[1], NULL );

  if( h == NULL) {
    fprintf( stderr, "Unable to Open file: %s\n", argv[1] );
    exit(2);
  }

  find_hdrval( h, "header_size", &header_val );
  header_size = *(int *)header_val.value;
  find_hdrval( h, "num_lags", &header_val );
  num_lags = *(int *)header_val.value;
  find_hdrval( h, "nifs", &header_val );     /*-- Number ifs? --*/
  nifs = *(int *)header_val.value;
  find_hdrval( h, "lagformat", &header_val );
  lagformat = *(int *)header_val.value;
  find_hdrval( h, "src_ra", &header_val );
  ra = *(double *)header_val.value;
  find_hdrval( h, "src_dec", &header_val );
  dec = *(double *)header_val.value;
  find_hdrval( h, "isalfa", &isalfa_val );
  isalfa = *(int *)isalfa_val.value;
  find_hdrval( h, "frontend", &header_val);
  strcpy(frontend, header_val.value);

  //Old
  //wapp = atoi( &((argv[1])[4]) );
  wapp = find_wapp(argv[1]);

  find_hdrval( h, "alfa_raj", &header_val );
  beam_raj = (double *)header_val.value;
  if(beam_raj == NULL)
    printf("Beam coordinates not found\n");
  else
    {
      find_hdrval( h, "alfa_decj", &header_val );
      beam_decj = (double *)header_val.value;
      find_hdrval( h, "alfa_az", &header_val );
      beam_az = (double *)header_val.value;
      find_hdrval( h, "alfa_za", &header_val );
      beam_za = (double *)header_val.value;
    }

  if( lagformat == 0 )
    datasz = sizeof(unsigned short);
  else
    datasz = sizeof(unsigned int);

  datalen = nifs * num_lags * datasz * 2;
  data = ( unsigned char *) malloc( 4*datalen );

  if( isalfa != 1 ) {
    fprintf( stderr, "According to header, not a dual pixel ALFA file %s\n isalfa = %d  frontend = %s\n", argv[1],isalfa,frontend );
    fprintf( stderr, "Splitting anyway...\n");
    //exit(2);
  }
  new_filename( argv[1], h, name );

  for( i=0; i<2; i++ ) {
    if( !(fd[i] = fopen( name[i], "w+" )) < 0 ) {
      fprintf( stderr, "unable to create %s\n", name[i] );
      continue;
    }
    printf( "opening %s\n", name[i] );
    
    if( fwrite( h->buf, h->offset, 1, fd[i] ) != 1 ) {
      fprintf( stderr, "unable to write ascii %s\n", name[i] );
      continue;
    }

    /* set isalfa to include alfa pixel number */
    //Old--gives a very large number
    //*(int *)isalfa_val.value = 2 |((wapp*2+i)<<4); 
    *(int *)isalfa_val.value = (wapp-1)*2 + i;
    /* if( (wapp*2 + i) != 0 ) */

    if(beam_raj == NULL)
      {	
	fix_positions( ra, dec, h, (wapp-1)*2+i );
      }
    else
      {
	ra = deg2dms(beam_raj[(wapp-1)*2+i]);
	find_hdrval( h, "src_ra", &header_val );
	*(double *)header_val.value = ra;
	
	dec = deg2dms(beam_decj[(wapp-1)*2+i]);
	find_hdrval( h, "src_dec", &header_val );
	*(double *)header_val.value = dec;
	/*
	find_hdrval( h, "start_az", &header_val );
	*(double *)header_val.value = beam_az[(wapp-1)*2+i];
	
	find_hdrval( h, "start_za", &header_val );
	*(double *)header_val.value = beam_za[(wapp-1)*2+i];
	*/
      }
    

    if( fwrite( h->header, h->headlen, 1, fd[i] ) != 1 ) {
      fprintf( stderr, "unable to write header %s\n", name[i] );
      continue;
    }
  }

  off = h->headlen + h->offset;
  lseek( h->fd, off, SEEK_SET );

  rdatalen = datalen;
  wdatalen = rdatalen/2;
  while(1) {
    if( ( n = read( h->fd, data, rdatalen )) == rdatalen ) {
      if( fwrite( data, wdatalen, 1, fd[0] ) != 1 ) {
        fprintf( stderr, "unable to write data %s\n", name[i] );
        break;
      }
      if( fwrite( &data[wdatalen], wdatalen, 1, fd[1] )!=1) {
        fprintf( stderr, "unable to write data %s\n", name[i] );
        break;
      }
    } else
      break;
    off += datalen;
  }
  fclose(fd[0]);
  fclose(fd[1]);
}

//Find wapp number regardless of the filename (which shouldn't contain
//any other w's)
int find_wapp(char filename[])
{
  int i = 0;
  while(filename[i] != 'w')
    i++;
  
  return atoi(&(filename[i+4]));
  
}


new_filename( old, h, name )
char *old;
struct HEADERP *h;
char *name[];
{
  char project[25];
  char obs_date[25];
  char start_time[25];
  char src_name[25];
  int mjd;
  int secs;
  int ix;
  char *source;
  int wapp;
  struct HEADERVAL header_val;

  find_hdrval( h, "project_id", &header_val );
  memcpy( project,  (char *)header_val.value, 24);
  project[24] = 0;

  find_hdrval( h, "obs_date", &header_val );
  memcpy( obs_date, (char *)header_val.value, 24);
  obs_date[24] = 0;
  mjd = obs2mjd(obs_date);

  find_hdrval( h, "start_time", &header_val );
  memcpy( start_time, (char *)header_val.value, 24);
  start_time[24] = 0;
  
  ix = find_index(old);
  wapp = find_wappno( old );

  find_hdrval( h, "src_name", &header_val );
  memcpy( src_name, (char *)header_val.value, 24);
  src_name[24] = 0;
  secs = ast_seconds(start_time);

  sprintf(name[0], "%s_%d_%05d_%04d_%s_%d.wapp",
	  project, (int)mjd, secs, ix, src_name, wapp*2 );

  sprintf(name[1], "%s_%d_%05d_%04d_%s_%d.wapp",
	  project, (int)mjd, secs, ix, src_name, wapp*2+1 );

}


int find_index(s)
char *s;
{
  s = strrchr( s, '.' );
  if( s )
    return atoi(&s[1]);
  else
    return 0;
}

int find_wappno(s)
char *s;
{
  int ret;

  ret = 0;
  if((s = strchr( s, '.' )))
    if((s = strstr( s, "wapp" )))
      ret = atoi(&s[4])-1;

  return ret;
}

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

  //fprintf(stderr,"cal2mjd call in obs2mjd: %f\n", cal2mjd( year, mon, day ));

  return( cal2mjd( year, mon, day ));
}

/*
   take the center position normally in src_ra, src_dec 
   and fix the positions with the ALFA offsets, from Desh
*/

void fix_positions( ra, dec, h, beam )
double ra;
double dec;
struct HEADERP *h;
int beam;
{
  double saz, sza, ast, nra, ndec;
  double lst, current_time, mjd, offx, offy, rot_angle;
  struct HEADERVAL header_val;
  int secs;
  char start_time[80];
  char obs_date[80];

/*
  if( beam <= 0 || beam >= 7 )
    return;
*/

  find_hdrval( h, "start_time", &header_val );
  memcpy( start_time, (char *)header_val.value, 24);
  secs = ast_seconds(start_time);

  find_hdrval( h, "obs_date", &header_val );
  memcpy( obs_date, (char *)header_val.value, 24);
  obs_date[24] = 0;
  mjd = obs2mjd(obs_date);
  mjd = mjd - obs2mjd("20000101");
  current_time = 2000.0 + (mjd + (secs/3600.0 - 4.0)/24.0)/365.25; 
   /* ast to ut */

  find_hdrval( h, "start_lst", &header_val );
  lst = *(double *)header_val.value;

  find_hdrval( h, "alfa_angle", &header_val );
  if( header_val.value )
    rot_angle = *(double *)header_val.value;
  else
    rot_angle = 0.0;

  ra = dms2deg(ra);
  dec = dms2deg(dec);

  alfa_position( ra, dec, lst, current_time, rot_angle, 
    0.0, 0.0, beam, &nra, &ndec );
  
  ra = deg2dms(nra);
  find_hdrval( h, "src_ra", &header_val );
  *(double *)header_val.value = ra;

  dec = deg2dms(ndec);
  find_hdrval( h, "src_dec", &header_val );
  *(double *)header_val.value = dec;
}


/*
*2.0*M_PI/24.0;
*M_PI/180.0;
  
  sla_preces("FK5",2000.0,current_time, &ra, &dec );
  ra = ra * 24.0 / 2.0 / M_PI;
  dec = dec * 180.0 / M_PI;

  compute_azza( ra, dec, lst, &saz, &sza );
  offx = find_alfaoff( beam, 0 );
  offy = find_alfaoff( beam, 1 );
  saz += rotx(offx, offy, rot_angle);
  sza += roty(offx, offy, rot_angle);
  compute_radec( saz, sza, lst, &ra, &dec );

  ra = ra*2.0*M_PI/24.0;
  dec = dec*M_PI/180.0;
  sla_preces("FK5",current_time, 2000.0, &ra, &dec );
  ra = ra * 24.0 / 2.0 / M_PI;
  dec = dec * 180.0 / M_PI;
*/

/* convert to  hhmmss.s */

double deg2dms( t )
double t;
{
  char buf[24];
  int hh, mm;
  double ss;
  int neg;

  neg =  t < 0;
  t = fabs(t);

  hh = (int)t;
  t = t - (double)hh;
  t = t * 60.0;
  mm = (int)t;
  t = t - (double)mm;

  ss = t*60.0;

  if (ss >= 60.0) {
    ss -= 60.0;
    mm++;
  }

  if( mm >= 60 ) {
    mm -= 60;
    hh++;
  }

  return( (neg ? -1.0 : 1.0) * (hh * 10000.0 + mm * 100.0 + ss) );
}

double dms2deg(in)
double in;
{
  int hh, mm, i, neg;
  double ss, ret;
  double si;

  si = fabs(in);
  neg = (in < 0.0 );

  hh = (int)(si/10000.0);
  si = si - hh * 10000.0;
  mm = (int)(si/100.0);
  ss = si - mm * 100.0;

  ret = abs(hh) + mm/60.0 + ss/3600.0;

  if( neg )
    ret = -ret;

  return( ret );
}


int ast_seconds(s)
char *s;
{
  int hh, mm, ss, ret;

  hh = atoi(s);
  mm = atoi( &s[3] );
  ss = atoi( &s[6] );

  ret = hh*3600 + mm * 60 + ss;
  
  return(ret);
}

