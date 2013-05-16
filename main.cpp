/****************************************************************************\
*	   Signal Server 1.3.4: Server optimised SPLAT! by Alex Farrant      *
******************************************************************************
*	SPLAT! Project started in 1997 by John A. Magliacane, KD2BD 	     *
*					                                     *
******************************************************************************
*         Please consult the SPLAT! documentation for a complete list of     *
*	     individuals who have contributed to this project. 		     *
******************************************************************************
*                                                                            *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.								     *
* 									     *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.							     *
*									     *
******************************************************************************
* g++ -Wall -O3 -s -lm -fomit-frame-pointer itm.cpp main.cpp -o ss           * 
\****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#define GAMMA 2.5
#define MAXPAGES 36
#define ARRAYSIZE 129650
#define IPPD 3600

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI 6.283185307179586
#endif

#ifndef HALFPI
#define HALFPI 1.570796326794896
#endif

#define DEG2RAD 1.74532925199e-02
#define EARTHRADIUS 20902230.97
#define	METERS_PER_MILE 1609.344
#define METERS_PER_FOOT 0.3048
#define	KM_PER_MILE 1.609344
#define FOUR_THIRDS 1.3333333333333

char 	string[255], sdf_path[255], udt_file[255], opened=0, gpsav=0, ss_name[16],
	ss_version[6], dashes[80];

double	earthradius, max_range=0.0, forced_erp, dpp, ppd,
	fzone_clearance=0.6, forced_freq, clutter, lat,lon,txh,tercon,terdic,
        north,east,south,west;

int	min_north=90, max_north=-90, min_west=360, max_west=-1, ippd, mpi,
	max_elevation=-32768, min_elevation=32768, bzerror, contour_threshold,
        pred,pblue,pgreen,ter,multiplier=256,debug=0,loops=64,jgets=0, MAXRAD;

unsigned char got_elevation_pattern, got_azimuth_pattern, metric=0, dbm=0;



struct site {	double lat;
		double lon;
		float alt;
		char name[50];
		char filename[255];
	    } 	site;

struct path {	double lat[ARRAYSIZE];
		double lon[ARRAYSIZE];
		double elevation[ARRAYSIZE];
		double distance[ARRAYSIZE];
		int length;
	    }	path;

struct dem {	int min_north;
		int max_north;
		int min_west;
		int max_west;
		int max_el;
		int min_el;
		short data[IPPD][IPPD];
		unsigned char mask[IPPD][IPPD];
		unsigned char signal[IPPD][IPPD];
           }	dem[MAXPAGES];

struct LR {	double eps_dielect;
		double sgm_conductivity;
		double eno_ns_surfref;
		double frq_mhz;
		double conf;
		double rel;
		double erp;
		int radio_climate;
		int pol;
		float antenna_pattern[361][1001];
          }	LR;

struct region { unsigned char color[128][3];
		int level[128];
		int levels;
	      }	region;

double elev[ARRAYSIZE+10];


void point_to_point(double elev[], double tht_m, double rht_m,
	  double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
	  double frq_mhz, int radio_climate, int pol, double conf,
	  double rel, double &dbloss, char *strmode, int &errnum);

double arccos(double x, double y)
{
	/* This function implements the arc cosine function,
	   returning a value between 0 and TWOPI. */

	double result=0.0;

	if (y>0.0)
		result=acos(x/y);

	if (y<0.0)
		result=PI+acos(x/y);

	return result;
}


int ReduceAngle(double angle)
{
	/* This function normalizes the argument to
	   an integer angle between 0 and 180 degrees */

	double temp;

	temp=acos(cos(angle*DEG2RAD));

	return (int)rint(temp/DEG2RAD);
}

double LonDiff(double lon1, double lon2)
{
	/* This function returns the short path longitudinal
	   difference between longitude1 and longitude2
	   as an angle between -180.0 and +180.0 degrees.
	   If lon1 is west of lon2, the result is positive.
	   If lon1 is east of lon2, the result is negative. */

	double diff;

	diff=lon1-lon2;

	if (diff<=-180.0)
		diff+=360.0;

	if (diff>=180.0)
		diff-=360.0;

	return diff;
}

char *dec2dms(double decimal)
{
	/* Converts decimal degrees to degrees, minutes, seconds,
	   (DMS) and returns the result as a character string. */

	char	sign;
	int	degrees, minutes, seconds;
	double	a, b, c, d;

	if (decimal<0.0)
	{
		decimal=-decimal;
		sign=-1;
	}

	else
		sign=1;

	a=floor(decimal);
	b=60.0*(decimal-a);
	c=floor(b);
	d=60.0*(b-c);

	degrees=(int)a;
	minutes=(int)c;
	seconds=(int)d;

	if (seconds<0)
		seconds=0;

	if (seconds>59)
		seconds=59;

	string[0]=0;
	snprintf(string,250,"%d%c %d\' %d\"", degrees*sign, 176, minutes, seconds);
	return (string);
}

int PutMask(double lat, double lon, int value)
{
	/* Lines, text, markings, and coverage areas are stored in a
	   mask that is combined with topology data when topographic
	   maps are generated by ss.  This function sets and resets
	   bits in the mask based on the latitude and longitude of the
	   area pointed to. */

	int	x, y, indx;
	char	found;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
	{
		dem[indx].mask[x][y]=value;
		return ((int)dem[indx].mask[x][y]);
	}

	else
		return -1;
}

int OrMask(double lat, double lon, int value)
{
	/* Lines, text, markings, and coverage areas are stored in a
	   mask that is combined with topology data when topographic
	   maps are generated by ss.  This function sets bits in
	   the mask based on the latitude and longitude of the area
	   pointed to. */

	int	x, y, indx;
	char	found;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
	{
		dem[indx].mask[x][y]|=value;
		return ((int)dem[indx].mask[x][y]);
	}

	else
		return -1;
}

int GetMask(double lat, double lon)
{
	/* This function returns the mask bits based on the latitude
	   and longitude given. */

	return (OrMask(lat,lon,0));
}

int PutSignal(double lat, double lon, unsigned char signal)
{
	/* This function writes a signal level (0-255)
	   at the specified location for later recall. */

	int	x, y, indx;
	char	found;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
	{
		dem[indx].signal[x][y]=signal;
		return (dem[indx].signal[x][y]);
	}

	else
		return 0;
}

unsigned char GetSignal(double lat, double lon)
{
	/* This function reads the signal level (0-255) at the
	   specified location that was previously written by the
	   complimentary PutSignal() function. */

	int	x, y, indx;
	char	found;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
		return (dem[indx].signal[x][y]);
	else
		return 0;
}

double GetElevation(struct site location)
{
	/* This function returns the elevation (in feet) of any location
	   represented by the digital elevation model data in memory.
	   Function returns -5000.0 for locations not found in memory. */

	char	found;
	int	x, y, indx;
	double	elevation;

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(location.lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,location.lon)));

		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
		elevation=3.28084*dem[indx].data[x][y];
	else
		elevation=-5000.0;

	return elevation;
}

int AddElevation(double lat, double lon, double height)
{
	/* This function adds a user-defined terrain feature
	   (in meters AGL) to the digital elevation model data
	   in memory.  Does nothing and returns 0 for locations
	   not found in memory. */

	char	found;
	int	x, y, indx;
	

	for (indx=0, found=0; indx<MAXPAGES && found==0;)
	{
		x=(int)rint(ppd*(lat-dem[indx].min_north));
		y=mpi-(int)rint(ppd*(LonDiff(dem[indx].max_west,lon)));
		
		if (x>=0 && x<=mpi && y>=0 && y<=mpi)
			found=1;
		else
			indx++;
	}

	if (found)
		dem[indx].data[x][y]+=(short)rint(height);

	return found;
}

double Distance(struct site site1, struct site site2)
{
	/* This function returns the great circle distance
	   in miles between any two site locations. */

	double	lat1, lon1, lat2, lon2, distance;

	lat1=site1.lat*DEG2RAD;
	lon1=site1.lon*DEG2RAD;
	lat2=site2.lat*DEG2RAD;
	lon2=site2.lon*DEG2RAD;

	distance=3959.0*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos((lon1)-(lon2)));

	return distance;
}

double Azimuth(struct site source, struct site destination)
{
	/* This function returns the azimuth (in degrees) to the
	   destination as seen from the location of the source. */

	double	dest_lat, dest_lon, src_lat, src_lon,
		beta, azimuth, diff, num, den, fraction;

	dest_lat=destination.lat*DEG2RAD;
	dest_lon=destination.lon*DEG2RAD;

	src_lat=source.lat*DEG2RAD;
	src_lon=source.lon*DEG2RAD;

	/* Calculate Surface Distance */

	beta=acos(sin(src_lat)*sin(dest_lat)+cos(src_lat)*cos(dest_lat)*cos(src_lon-dest_lon));

	/* Calculate Azimuth */

	num=sin(dest_lat)-(sin(src_lat)*cos(beta));
	den=cos(src_lat)*sin(beta);
	fraction=num/den;

	/* Trap potential problems in acos() due to rounding */

	if (fraction>=1.0)
		fraction=1.0;

	if (fraction<=-1.0)
		fraction=-1.0;

	/* Calculate azimuth */

	azimuth=acos(fraction);

	/* Reference it to True North */

	diff=dest_lon-src_lon;

	if (diff<=-PI)
		diff+=TWOPI;

	if (diff>=PI)
		diff-=TWOPI;

	if (diff>0.0)
		azimuth=TWOPI-azimuth;

	return (azimuth/DEG2RAD);
}

double ElevationAngle(struct site source, struct site destination)
{
	/* This function returns the angle of elevation (in degrees)
	   of the destination as seen from the source location.
	   A positive result represents an angle of elevation (uptilt),
	   while a negative result represents an angle of depression
	   (downtilt), as referenced to a normal to the center of
	   the earth. */

	register double a, b, dx;

	a=GetElevation(destination)+destination.alt+earthradius;
	b=GetElevation(source)+source.alt+earthradius;

 	dx=5280.0*Distance(source,destination);

	/* Apply the Law of Cosines */

	return ((180.0*(acos(((b*b)+(dx*dx)-(a*a))/(2.0*b*dx)))/PI)-90.0);
}

void ReadPath(struct site source, struct site destination)
{
	/* This function generates a sequence of latitude and
	   longitude positions between source and destination
	   locations along a great circle path, and stores
	   elevation and distance information for points
	   along that path in the "path" structure. */

	int	c;
	double	azimuth, distance, lat1, lon1, beta, den, num,
		lat2, lon2, total_distance, dx, dy, path_length,
		miles_per_sample, samples_per_radian=68755.0;
	struct	site tempsite;
	
	lat1=source.lat*DEG2RAD;
	lon1=source.lon*DEG2RAD;

	lat2=destination.lat*DEG2RAD;
	lon2=destination.lon*DEG2RAD;

        samples_per_radian=ppd*57.295833;
        
       	azimuth=Azimuth(source,destination)*DEG2RAD;

	total_distance=Distance(source,destination);

	if (total_distance>(30.0/ppd))		
	{
		dx=samples_per_radian*acos(cos(lon1-lon2));
		dy=samples_per_radian*acos(cos(lat1-lat2));

		path_length=sqrt((dx*dx)+(dy*dy));		

		miles_per_sample=total_distance/path_length;	
	}

	else
	{
		c=0;
		dx=0.0;
		dy=0.0;
		path_length=0.0;
		miles_per_sample=0.0;
		total_distance=0.0;

		lat1=lat1/DEG2RAD;
		lon1=lon1/DEG2RAD;

		path.lat[c]=lat1;
		path.lon[c]=lon1;
		path.elevation[c]=GetElevation(source);
		path.distance[c]=0.0;
	}

	for (distance=0.0, c=0; (total_distance!=0.0 && distance<=total_distance && c<ARRAYSIZE); c++, distance=miles_per_sample*(double)c)
	{
		beta=distance/3959.0;
		lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
		num=cos(beta)-(sin(lat1)*sin(lat2));
		den=cos(lat1)*cos(lat2);

		if (azimuth==0.0 && (beta>HALFPI-lat1))
			lon2=lon1+PI;

		else if (azimuth==HALFPI && (beta>HALFPI+lat1))
				lon2=lon1+PI;

		else if (fabs(num/den)>1.0)
				lon2=lon1;

		else
		{
			if ((PI-azimuth)>=0.0)
				lon2=lon1-arccos(num,den);
			else
				lon2=lon1+arccos(num,den);
		}

		while (lon2<0.0)
			lon2+=TWOPI;

		while (lon2>TWOPI)
			lon2-=TWOPI;

		lat2=lat2/DEG2RAD;
		lon2=lon2/DEG2RAD;

		path.lat[c]=lat2;
		path.lon[c]=lon2;
		tempsite.lat=lat2;
		tempsite.lon=lon2;
		path.elevation[c]=GetElevation(tempsite);
		path.distance[c]=distance;
	}

	/* Make sure exact destination point is recorded at path.length-1 */

	if (c<ARRAYSIZE)
	{
		path.lat[c]=destination.lat;
		path.lon[c]=destination.lon;
		path.elevation[c]=GetElevation(destination);
		path.distance[c]=total_distance;
		c++;
	}

	if (c<ARRAYSIZE)
		path.length=c;
	else
		path.length=ARRAYSIZE-1;
}

double ElevationAngle2(struct site source, struct site destination, double er)
{
	/* This function returns the angle of elevation (in degrees)
	   of the destination as seen from the source location, UNLESS
	   the path between the sites is obstructed, in which case, the
	   elevation angle to the first obstruction is returned instead.
	   "er" represents the earth radius. */

	int	x;
	char	block=0;
	double	source_alt, destination_alt, cos_xmtr_angle,
		cos_test_angle, test_alt, elevation, distance,
		source_alt2, first_obstruction_angle=0.0;
	struct	path temp;

	temp=path;

	ReadPath(source,destination);

	distance=5280.0*Distance(source,destination);
	source_alt=er+source.alt+GetElevation(source);
	destination_alt=er+destination.alt+GetElevation(destination);
	source_alt2=source_alt*source_alt;

	/* Calculate the cosine of the elevation angle of the
	   destination (receiver) as seen by the source (transmitter). */

	cos_xmtr_angle=((source_alt2)+(distance*distance)-(destination_alt*destination_alt))/(2.0*source_alt*distance);

	/* Test all points in between source and destination locations to
	   see if the angle to a topographic feature generates a higher
	   elevation angle than that produced by the destination.  Begin
	   at the source since we're interested in identifying the FIRST
	   obstruction along the path between source and destination. */

	for (x=2, block=0; x<path.length && block==0; x++)
	{
		distance=5280.0*path.distance[x];

		test_alt=earthradius+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

		cos_test_angle=((source_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*source_alt*distance);

		/* Compare these two angles to determine if
		   an obstruction exists.  Since we're comparing
		   the cosines of these angles rather than
		   the angles themselves, the sense of the
		   following "if" statement is reversed from
		   what it would be if the angles themselves
		   were compared. */

		if (cos_xmtr_angle>=cos_test_angle)
		{
			block=1;
			first_obstruction_angle=((acos(cos_test_angle))/DEG2RAD)-90.0;
		}
	}

	if (block)
		elevation=first_obstruction_angle;

	else
		elevation=((acos(cos_xmtr_angle))/DEG2RAD)-90.0;

	path=temp;

	return elevation;
}

double AverageTerrain(struct site source, double azimuthx, double start_distance, double end_distance)
{
	/* This function returns the average terrain calculated in
	   the direction of "azimuth" (degrees) between "start_distance"
	   and "end_distance" (miles) from the source location.  If
	   the terrain is all water (non-critical error), -5000.0 is
	   returned.  If not enough SDF data has been loaded into
	   memory to complete the survey (critical error), then
	   -9999.0 is returned. */

	int	c, samples, endpoint;
	double	beta, lat1, lon1, lat2, lon2, num, den, azimuth, terrain=0.0;
	struct	site destination;

	lat1=source.lat*DEG2RAD;
	lon1=source.lon*DEG2RAD;

	/* Generate a path of elevations between the source
	   location and the remote location provided. */

	beta=end_distance/3959.0;

	azimuth=DEG2RAD*azimuthx;

	lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
	num=cos(beta)-(sin(lat1)*sin(lat2));
	den=cos(lat1)*cos(lat2);

	if (azimuth==0.0 && (beta>HALFPI-lat1))
		lon2=lon1+PI;

	else if (azimuth==HALFPI && (beta>HALFPI+lat1))
			lon2=lon1+PI;

	else if (fabs(num/den)>1.0)
			lon2=lon1;

	else
	{
		if ((PI-azimuth)>=0.0)
			lon2=lon1-arccos(num,den);
		else
			lon2=lon1+arccos(num,den);
	}

	while (lon2<0.0)
		lon2+=TWOPI;

	while (lon2>TWOPI)
		lon2-=TWOPI;

	lat2=lat2/DEG2RAD;
	lon2=lon2/DEG2RAD;

	destination.lat=lat2;
	destination.lon=lon2;

	/* If SDF data is missing for the endpoint of
	   the radial, then the average terrain cannot
	   be accurately calculated.  Return -9999.0 */

	if (GetElevation(destination)<-4999.0)
		return (-9999.0);
	else
	{
		ReadPath(source,destination);

		endpoint=path.length;

		/* Shrink the length of the radial if the
		   outermost portion is not over U.S. land. */

		for (c=endpoint-1; c>=0 && path.elevation[c]==0.0; c--);

		endpoint=c+1;

		for (c=0, samples=0; c<endpoint; c++)
		{
			if (path.distance[c]>=start_distance)
			{
				terrain+=(path.elevation[c]==0.0?path.elevation[c]:path.elevation[c]+clutter);
				samples++;
			}
		}

		if (samples==0)
			terrain=-5000.0;  /* No land */
		else
			terrain=(terrain/(double)samples);

		return terrain;
	}
}

double haat(struct site antenna)
{
	/* This function returns the antenna's Height Above Average
	   Terrain (HAAT) based on FCC Part 73.313(d).  If a critical
	   error occurs, such as a lack of SDF data to complete the
	   survey, -5000.0 is returned. */

	int	azi, c;
	char	error=0;
	double	terrain, avg_terrain, haat, sum=0.0;

	/* Calculate the average terrain between 2 and 10 miles
	   from the antenna site at azimuths of 0, 45, 90, 135,
	   180, 225, 270, and 315 degrees. */

	for (c=0, azi=0; azi<=315 && error==0; azi+=45)
	{
		terrain=AverageTerrain(antenna, (double)azi, 2.0, 10.0);

		if (terrain<-9998.0)  /* SDF data is missing */
			error=1;

		if (terrain>-4999.0)  /* It's land, not water */
		{
			sum+=terrain;  /* Sum of averages */
			c++;
		}
	}

	if (error)
		return -5000.0;
	else
	{
		avg_terrain=(sum/(double)c);
		haat=(antenna.alt+GetElevation(antenna))-avg_terrain;
		return haat;
	}
}
double ReadBearing(char *input)
{
	/* This function takes numeric input in the form of a character
	   string, and returns an equivalent bearing in degrees as a
	   decimal number (double).  The input may either be expressed
	   in decimal format (40.139722) or degree, minute, second
	   format (40 08 23).  This function also safely handles
	   extra spaces found either leading, trailing, or
	   embedded within the numbers expressed in the
	   input string.  Decimal seconds are permitted. */

	double	seconds, bearing=0.0;
	char	string[20];
	int	a, b, length, degrees, minutes;

	/* Copy "input" to "string", and ignore any extra
	   spaces that might be present in the process. */

	string[0]=0;
	length=strlen(input);

	for (a=0, b=0; a<length && a<18; a++)
	{
		if ((input[a]!=32 && input[a]!='\n') || (input[a]==32 && input[a+1]!=32 && input[a+1]!='\n' && b!=0))
		{
			string[b]=input[a];
			b++;
		}
	}

	string[b]=0;

	/* Count number of spaces in the clean string. */

	length=strlen(string);

	for (a=0, b=0; a<length; a++)
		if (string[a]==32)
			b++;

	if (b==0)  /* Decimal Format (40.139722) */
		sscanf(string,"%lf",&bearing);

	if (b==2)  /* Degree, Minute, Second Format (40 08 23.xx) */
	{
		sscanf(string,"%d %d %lf",&degrees, &minutes, &seconds);

		bearing=fabs((double)degrees);
		bearing+=fabs(((double)minutes)/60.0);
		bearing+=fabs(seconds/3600.0);

		if ((degrees<0) || (minutes<0) || (seconds<0.0))
			bearing=-bearing;
	}

	/* Anything else returns a 0.0 */

	if (bearing>360.0 || bearing<-360.0)
		bearing=0.0;

	return bearing;
}

void LoadPAT(char *filename)
{
	/* This function reads and processes antenna pattern (.az
	   and .el) files that correspond in name to previously
	   loaded ss .lrp files.  */

	int	a, b, w, x, y, z, last_index, next_index, span;
	char	string[255], azfile[255], elfile[255], *pointer=NULL, *s=NULL;
	float	az, xx, elevation, amplitude, rotation, valid1, valid2,
		delta, azimuth[361], azimuth_pattern[361], el_pattern[10001],
		elevation_pattern[361][1001], slant_angle[361], tilt,
		mechanical_tilt=0.0, tilt_azimuth, tilt_increment, sum;
	FILE	*fd=NULL;
	unsigned char read_count[10001];

	for (x=0; filename[x]!='.' && filename[x]!=0 && x<250; x++)
	{
		azfile[x]=filename[x];
		elfile[x]=filename[x];
	}

	azfile[x]='.';
	azfile[x+1]='a';
	azfile[x+2]='z';
	azfile[x+3]=0;

	elfile[x]='.';
	elfile[x+1]='e';
	elfile[x+2]='l';
	elfile[x+3]=0;

	rotation=0.0;

	got_azimuth_pattern=0;
	got_elevation_pattern=0;

	/* Load .az antenna pattern file */

	fd=fopen(azfile,"r");

	if (fd!=NULL)
	{
		/* Clear azimuth pattern array */

		for (x=0; x<=360; x++)
		{
			azimuth[x]=0.0;
			read_count[x]=0;
		}


		/* Read azimuth pattern rotation
		   in degrees measured clockwise
		   from true North. */

		s=fgets(string,254,fd);
		pointer=strchr(string,';');

		if (pointer!=NULL)
			*pointer=0;

		sscanf(string,"%f",&rotation);


		/* Read azimuth (degrees) and corresponding
		   normalized field radiation pattern amplitude
		   (0.0 to 1.0) until EOF is reached. */

		s=fgets(string,254,fd);
		pointer=strchr(string,';');

		if (pointer!=NULL)
			*pointer=0;

		sscanf(string,"%f %f",&az, &amplitude);

		do
		{
			x=(int)rintf(az);

			if (x>=0 && x<=360 && fd!=NULL)
			{
				azimuth[x]+=amplitude;
				read_count[x]++;
			}

			s=fgets(string,254,fd);
			pointer=strchr(string,';');

			if (pointer!=NULL)
				*pointer=0;

			sscanf(string,"%f %f",&az, &amplitude);

		} while (feof(fd)==0);

		fclose(fd);


		/* Handle 0=360 degree ambiguity */

		if ((read_count[0]==0) && (read_count[360]!=0))
		{
			read_count[0]=read_count[360];
			azimuth[0]=azimuth[360];
		}

		if ((read_count[0]!=0) && (read_count[360]==0))
		{
			read_count[360]=read_count[0];
			azimuth[360]=azimuth[0];
		}

		/* Average pattern values in case more than
		    one was read for each degree of azimuth. */

		for (x=0; x<=360; x++)
		{
			if (read_count[x]>1)
				azimuth[x]/=(float)read_count[x];
		}

		/* Interpolate missing azimuths
		   to completely fill the array */

		last_index=-1;
		next_index=-1;

		for (x=0; x<=360; x++)
		{
			if (read_count[x]!=0)
			{
				if (last_index==-1)
					last_index=x;
				else
					next_index=x;
			}

			if (last_index!=-1 && next_index!=-1)
			{
				valid1=azimuth[last_index];
				valid2=azimuth[next_index];

				span=next_index-last_index;
				delta=(valid2-valid1)/(float)span;

				for (y=last_index+1; y<next_index; y++)
					azimuth[y]=azimuth[y-1]+delta;

				last_index=y;
				next_index=-1;
			}
		}

		/* Perform azimuth pattern rotation
		   and load azimuth_pattern[361] with
		   azimuth pattern data in its final form. */

		for (x=0; x<360; x++)
		{
			y=x+(int)rintf(rotation);

			if (y>=360)
				y-=360;

			azimuth_pattern[y]=azimuth[x];
		}

		azimuth_pattern[360]=azimuth_pattern[0];

		got_azimuth_pattern=255;
	}

	/* Read and process .el file */

	fd=fopen(elfile,"r");

	if (fd!=NULL)
	{
		for (x=0; x<=10000; x++)
		{
			el_pattern[x]=0.0;
			read_count[x]=0;
		}

		/* Read mechanical tilt (degrees) and
		   tilt azimuth in degrees measured
		   clockwise from true North. */

		s=fgets(string,254,fd);
		pointer=strchr(string,';');

		if (pointer!=NULL)
			*pointer=0;

		sscanf(string,"%f %f",&mechanical_tilt, &tilt_azimuth);

		/* Read elevation (degrees) and corresponding
		   normalized field radiation pattern amplitude
		   (0.0 to 1.0) until EOF is reached. */

		s=fgets(string,254,fd);
		pointer=strchr(string,';');

		if (pointer!=NULL)
			*pointer=0;

		sscanf(string,"%f %f", &elevation, &amplitude);

		while (feof(fd)==0)
		{
			/* Read in normalized radiated field values
			   for every 0.01 degrees of elevation between
			   -10.0 and +90.0 degrees */

			x=(int)rintf(100.0*(elevation+10.0));

			if (x>=0 && x<=10000)
			{
				el_pattern[x]+=amplitude;
				read_count[x]++;
			}

			s=fgets(string,254,fd);
			pointer=strchr(string,';');

			if (pointer!=NULL)
				*pointer=0;

			sscanf(string,"%f %f", &elevation, &amplitude);
		}

		fclose(fd);

		/* Average the field values in case more than
		   one was read for each 0.01 degrees of elevation. */

		for (x=0; x<=10000; x++)
		{
			if (read_count[x]>1)
				el_pattern[x]/=(float)read_count[x];
		}

		/* Interpolate between missing elevations (if
		   any) to completely fill the array and provide
		   radiated field values for every 0.01 degrees of
		   elevation. */

		last_index=-1;
		next_index=-1;

		for (x=0; x<=10000; x++)
		{
			if (read_count[x]!=0)
			{
				if (last_index==-1)
					last_index=x;
				else
					next_index=x;
			}

			if (last_index!=-1 && next_index!=-1)
			{
				valid1=el_pattern[last_index];
				valid2=el_pattern[next_index];

				span=next_index-last_index;
				delta=(valid2-valid1)/(float)span;

				for (y=last_index+1; y<next_index; y++)
					el_pattern[y]=el_pattern[y-1]+delta;

				last_index=y;
				next_index=-1;
			}
		}

		/* Fill slant_angle[] array with offset angles based
		   on the antenna's mechanical beam tilt (if any)
		   and tilt direction (azimuth). */

		if (mechanical_tilt==0.0)
		{
			for (x=0; x<=360; x++)
				slant_angle[x]=0.0;
		}

		else
		{
			tilt_increment=mechanical_tilt/90.0;

			for (x=0; x<=360; x++)
			{
				xx=(float)x;
				y=(int)rintf(tilt_azimuth+xx);

				while (y>=360)
					y-=360;

				while (y<0)
					y+=360;

				if (x<=180)
					slant_angle[y]=-(tilt_increment*(90.0-xx));

				if (x>180)
					slant_angle[y]=-(tilt_increment*(xx-270.0));
			}
		}

		slant_angle[360]=slant_angle[0];   /* 360 degree wrap-around */

		for (w=0; w<=360; w++)
		{
			tilt=slant_angle[w];

			/** Convert tilt angle to
			    an array index offset **/

			y=(int)rintf(100.0*tilt);

			/* Copy shifted el_pattern[10001] field
			   values into elevation_pattern[361][1001]
			   at the corresponding azimuth, downsampling
			   (averaging) along the way in chunks of 10. */

			for (x=y, z=0; z<=1000; x+=10, z++)
			{
				for (sum=0.0, a=0; a<10; a++)
				{
					b=a+x;

					if (b>=0 && b<=10000)
						sum+=el_pattern[b];
					if (b<0)
						sum+=el_pattern[0];
					if (b>10000)
						sum+=el_pattern[10000];
				}

				elevation_pattern[w][z]=sum/10.0;
			}
		}

		got_elevation_pattern=255;
	}

	for (x=0; x<=360; x++)
	{
		for (y=0; y<=1000; y++)
		{
			if (got_elevation_pattern)
				elevation=elevation_pattern[x][y];
			else
				elevation=1.0;

			if (got_azimuth_pattern)
				az=azimuth_pattern[x];
			else
				az=1.0;

			LR.antenna_pattern[x][y]=az*elevation;
		}
	}
}

int LoadSDF_SDF(char *name)
{
	/* This function reads uncompressed ss Data Files (.sdf)
	   containing digital elevation model data into memory.
	   Elevation data, maximum and minimum elevations, and
	   quadrangle limits are stored in the first available
	   dem[] structure. */

	int	x, y, data, indx, minlat, minlon, maxlat, maxlon,j;
	char	found, free_page=0, line[20], jline[20], sdf_file[255],
		path_plus_name[255], *s=NULL,*junk=NULL;
        
        
	FILE	*fd;

	for (x=0; name[x]!='.' && name[x]!=0 && x<250; x++)
		sdf_file[x]=name[x];

	sdf_file[x]=0;

	/* Parse filename for minimum latitude and longitude values */

	sscanf(sdf_file,"%d:%d:%d:%d",&minlat,&maxlat,&minlon,&maxlon);

	sdf_file[x]='.';
	sdf_file[x+1]='s';
	sdf_file[x+2]='d';
	sdf_file[x+3]='f';
	sdf_file[x+4]=0;

	/* Is it already in memory? */
	

	for (indx=0, found=0; indx<MAXPAGES && found==0; indx++)
	{
		if (minlat==dem[indx].min_north && minlon==dem[indx].min_west && maxlat==dem[indx].max_north && maxlon==dem[indx].max_west)
			found=1;
	}

	/* Is room available to load it? */

	if (found==0)
	{
		for (indx=0, free_page=0; indx<MAXPAGES && free_page==0; indx++)
			if (dem[indx].max_north==-90)
				free_page=1;
	}

	indx--;

	if (free_page && found==0 && indx>=0 && indx<MAXPAGES)
	{
		/* Search for SDF file in current working directory first */

		strncpy(path_plus_name,sdf_file,255);
		
			
		fd=fopen(path_plus_name,"rb");

		if (fd==NULL)
		{
			/* Next, try loading SDF file from path specified
			   in $HOME/.ss_path file or by -d argument */

			strncpy(path_plus_name,sdf_path,255);
			strncat(path_plus_name,sdf_file,255);
			
			

			fd=fopen(path_plus_name,"rb");
		}

		if (fd!=NULL)
		{
                    if (debug==1){
			fprintf(stdout,"Loading \"%s\" into page %d...",path_plus_name,indx+1);
			fflush(stdout);
                    }

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].max_west);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].min_north);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].min_west);

			s=fgets(line,19,fd);
			sscanf(line,"%d",&dem[indx].max_north);

                        /*
                         Here X lines of DEM will be read until IPPD is reached.
                         Each .sdf tile contains 1200x1200 = 1.44M 'points'
                         Each point is sampled for 1200 resolution!
                         
                        
                         */
			for (x=0; x<ippd; x++)
                        {
				for (y=0; y<ippd; y++)
				{
                            
                                      for (j=0; j<jgets; j++)
                                      {
                                        junk=fgets(jline,19,fd);   
                                        }
                                    
                                        s=fgets(line,19,fd);
                                                                        
					data=atoi(line); 
                                        
                                    	dem[indx].data[x][y]=data;
					dem[indx].signal[x][y]=0;
					dem[indx].mask[x][y]=0;

					if (data>dem[indx].max_el)
						dem[indx].max_el=data;

					if (data<dem[indx].min_el)
						dem[indx].min_el=data;
                                       
				} 
                                
                                if (ippd==600){ 
                                for (j=0; j<IPPD; j++)
                                      {
                                        junk=fgets(jline,19,fd);   
                                        }
                                }
                                if (ippd==300){ 
                                for (j=0; j<IPPD; j++)
                                      {
                                        junk=fgets(jline,19,fd); 
                                        junk=fgets(jline,19,fd);
                                        junk=fgets(jline,19,fd);
                                                                           
                                        }
                                }
                                
                        } 
                        
			fclose(fd);
                                        
			if (dem[indx].min_el<min_elevation)
				min_elevation=dem[indx].min_el;

			if (dem[indx].max_el>max_elevation)
				max_elevation=dem[indx].max_el;

			if (max_north==-90)
				max_north=dem[indx].max_north;

			else if (dem[indx].max_north>max_north)
					max_north=dem[indx].max_north;

			if (min_north==90)
				min_north=dem[indx].min_north;

			else if (dem[indx].min_north<min_north)
					min_north=dem[indx].min_north;

			if (max_west==-1)
				max_west=dem[indx].max_west;

			else
			{
				if (abs(dem[indx].max_west-max_west)<180)
				{
 					if (dem[indx].max_west>max_west)
						max_west=dem[indx].max_west;
				}

				else
				{
 					if (dem[indx].max_west<max_west)
						max_west=dem[indx].max_west;
				}
			}

			if (min_west==360)
				min_west=dem[indx].min_west;

			else
			{
				if (fabs(dem[indx].min_west-min_west)<180.0)
				{
 					if (dem[indx].min_west<min_west)
						min_west=dem[indx].min_west;
				}

				else
				{
 					if (dem[indx].min_west>min_west)
						min_west=dem[indx].min_west;
				}
			}

                       
			return 1;
		}

		else
			return -1;
	}

	else
		return 0;
}
char LoadSDF(char *name)
{
	/* This function loads the requested SDF file from the filesystem.
	   It first tries to invoke the LoadSDF_SDF() function to load an
	   uncompressed SDF file (since uncompressed files load slightly
	   faster).  If that attempt fails, then it tries to load a
	   compressed SDF file by invoking the LoadSDF_BZ() function.
	   If that fails, then we can assume that no elevation data
	   exists for the region requested, and that the region
	   requested must be entirely over water. */

	int	x, y, indx, minlat, minlon, maxlat, maxlon;
	char	found, free_page=0;
	int	return_value=-1;

	return_value=LoadSDF_SDF(name);


	/* If neither format can be found, then assume the area is water. */

	if (return_value==0 || return_value==-1)
	{
			
    			
	
		/* Parse SDF name for minimum latitude and longitude values */

		sscanf(name,"%d:%d:%d:%d",&minlat,&maxlat,&minlon,&maxlon);

		/* Is it already in memory? */

		for (indx=0, found=0; indx<MAXPAGES && found==0; indx++)
		{
			if (minlat==dem[indx].min_north && minlon==dem[indx].min_west && maxlat==dem[indx].max_north && maxlon==dem[indx].max_west)
				found=1;
		}

		/* Is room available to load it? */

		if (found==0)
		{
			for (indx=0, free_page=0; indx<MAXPAGES && free_page==0; indx++)
				if (dem[indx].max_north==-90)
					free_page=1;
		}

		indx--;

		if (free_page && found==0 && indx>=0 && indx<MAXPAGES)
		{
                    if(debug==1){
			fprintf(stdout,"Region  \"%s\" assumed as sea-level into page %d...",name,indx+1);
			fflush(stdout);
                    }

			dem[indx].max_west=maxlon;
			dem[indx].min_north=minlat;
			dem[indx].min_west=minlon;
			dem[indx].max_north=maxlat;

			/* Fill DEM with sea-level topography */

			for (x=0; x<ippd; x++)
				for (y=0; y<ippd; y++)
				{
		    			dem[indx].data[x][y]=0;
					dem[indx].signal[x][y]=0;
					dem[indx].mask[x][y]=0;

					if (dem[indx].min_el>0)
						dem[indx].min_el=0;
				}

			if (dem[indx].min_el<min_elevation)
				min_elevation=dem[indx].min_el;

			if (dem[indx].max_el>max_elevation)
				max_elevation=dem[indx].max_el;

			if (max_north==-90)
				max_north=dem[indx].max_north;

			else if (dem[indx].max_north>max_north)
					max_north=dem[indx].max_north;

			if (min_north==90)
				min_north=dem[indx].min_north;

			else if (dem[indx].min_north<min_north)
					min_north=dem[indx].min_north;

			if (max_west==-1)
				max_west=dem[indx].max_west;

			else
			{
				if (abs(dem[indx].max_west-max_west)<180)
				{
 					if (dem[indx].max_west>max_west)
						max_west=dem[indx].max_west;
				}

				else
				{
 					if (dem[indx].max_west<max_west)
						max_west=dem[indx].max_west;
				}
			}

			if (min_west==360)
				min_west=dem[indx].min_west;

			else
			{
				if (abs(dem[indx].min_west-min_west)<180)
				{
 					if (dem[indx].min_west<min_west)
						min_west=dem[indx].min_west;
				}

				else
				{
 					if (dem[indx].min_west>min_west)
						min_west=dem[indx].min_west;
				}
			}
                       
			return_value=1;
		}
	}

	return return_value;
}

void PlotPath(struct site source, struct site destination, char mask_value)
{
	/* This function analyzes the path between the source and
	   destination locations.  It determines which points along
	   the path have line-of-sight visibility to the source.
	   Points along with path having line-of-sight visibility
	   to the source at an AGL altitude equal to that of the
	   destination location are stored by setting bit 1 in the
	   mask[][] array, which are displayed in green when PPM
	   maps are later generated by ss. */

	char block;
	int x, y;
	register double cos_xmtr_angle, cos_test_angle, test_alt;
	double distance, rx_alt, tx_alt;

	ReadPath(source,destination);

	for (y=0; y<path.length; y++)
	{
		/* Test this point only if it hasn't been already
		   tested and found to be free of obstructions. */

		if ((GetMask(path.lat[y],path.lon[y])&mask_value)==0)
		{
			distance=5280.0*path.distance[y];
			tx_alt=earthradius+source.alt+path.elevation[0];
			rx_alt=earthradius+destination.alt+path.elevation[y];

			/* Calculate the cosine of the elevation of the
			   transmitter as seen at the temp rx point. */

			cos_xmtr_angle=((rx_alt*rx_alt)+(distance*distance)-(tx_alt*tx_alt))/(2.0*rx_alt*distance);

			for (x=y, block=0; x>=0 && block==0; x--)
			{
				distance=5280.0*(path.distance[y]-path.distance[x]);
				test_alt=earthradius+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

				cos_test_angle=((rx_alt*rx_alt)+(distance*distance)-(test_alt*test_alt))/(2.0*rx_alt*distance);

				/* Compare these two angles to determine if
				   an obstruction exists.  Since we're comparing
				   the cosines of these angles rather than
				   the angles themselves, the following "if"
				   statement is reversed from what it would
				   be if the actual angles were compared. */

				if (cos_xmtr_angle>=cos_test_angle)
					block=1;
			}

			if (block==0)
				OrMask(path.lat[y],path.lon[y],mask_value);
		}
	}
}

void PlotLRPath(struct site source, struct site destination, unsigned char mask_value, FILE *fd)
{
	/* This function plots the RF path loss between source and
	   destination points based on the Longley-Rice propagation
	   model, taking into account antenna pattern data, if available. */

	int	x, y, ifs, ofs, errnum;
	char	block=0, strmode[100];
	double	loss, azimuth, pattern=0.0,
		xmtr_alt, dest_alt, xmtr_alt2, dest_alt2,
		cos_rcvr_angle, cos_test_angle=0.0, test_alt,
		elevation=0.0, distance=0.0, four_thirds_earth,
		field_strength=0.0, rxp, dBm;
	struct	site temp;

	ReadPath(source,destination);

	four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;

	/* Copy elevations plus clutter along path into the elev[] array. */

	for (x=1; x<path.length-1; x++)
		elev[x+2]=(path.elevation[x]==0.0?path.elevation[x]*METERS_PER_FOOT:(clutter+path.elevation[x])*METERS_PER_FOOT);

	/* Copy ending points without clutter */

	elev[2]=path.elevation[0]*METERS_PER_FOOT;
	elev[path.length+1]=path.elevation[path.length-1]*METERS_PER_FOOT;

	/* Since the only energy the Longley-Rice model considers
	   reaching the destination is based on what is scattered
	   or deflected from the first obstruction along the path,
	   we first need to find the location and elevation angle
	   of that first obstruction (if it exists).  This is done
	   using a 4/3rds Earth radius to match the model used by
	   Longley-Rice.  This information is required for properly
	   integrating the antenna's elevation pattern into the
	   calculation for overall path loss. */

	for (y=2; (y<(path.length-1) && path.distance[y]<=max_range); y++)
	{
		/* Process this point only if it
		   has not already been processed. */
        
	
        if ((GetMask(path.lat[y],path.lon[y])&248)!=(mask_value<<3))
		{
			distance=5280.0*path.distance[y];
			xmtr_alt=four_thirds_earth+source.alt+path.elevation[0];
			dest_alt=four_thirds_earth+destination.alt+path.elevation[y];
			dest_alt2=dest_alt*dest_alt;
			xmtr_alt2=xmtr_alt*xmtr_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_rcvr_angle=((xmtr_alt2)+(distance*distance)-(dest_alt2))/(2.0*xmtr_alt*distance);

			if (cos_rcvr_angle>1.0)
				cos_rcvr_angle=1.0;

			if (cos_rcvr_angle<-1.0)
				cos_rcvr_angle=-1.0;

			if (got_elevation_pattern || fd!=NULL)
			{
				/* Determine the elevation angle to the first obstruction
				   along the path IF elevation pattern data is available
				   or an output (.ano) file has been designated. */

				for (x=2, block=0; (x<y && block==0); x++)
				{
					distance=5280.0*path.distance[x];

					test_alt=four_thirds_earth+(path.elevation[x]==0.0?path.elevation[x]:path.elevation[x]+clutter);

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle=((xmtr_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*xmtr_alt*distance);

					if (cos_test_angle>1.0)
						cos_test_angle=1.0;

					if (cos_test_angle<-1.0)
						cos_test_angle=-1.0;

					/* Compare these two angles to determine if
					   an obstruction exists.  Since we're comparing
					   the cosines of these angles rather than
					   the angles themselves, the sense of the
					   following "if" statement is reversed from
				  	   what it would be if the angles themselves
					   were compared. */

					if (cos_rcvr_angle>=cos_test_angle)
						block=1;
				}

				if (block)
					elevation=((acos(cos_test_angle))/DEG2RAD)-90.0;
				else
					elevation=((acos(cos_rcvr_angle))/DEG2RAD)-90.0;
			}

			/* Determine attenuation for each point along the
			   path using Longley-Rice's point_to_point mode
			   starting at y=2 (number_of_points = 1), the
			   shortest distance terrain can play a role in
			   path loss. */

			elev[0]=y-1;  /* (number of points - 1) */

			/* Distance between elevation samples */

			elev[1]=METERS_PER_MILE*(path.distance[y]-path.distance[y-1]);

			point_to_point(elev,source.alt*METERS_PER_FOOT,
	   		destination.alt*METERS_PER_FOOT, LR.eps_dielect,
			LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
			LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
			strmode, errnum);

			temp.lat=path.lat[y];
			temp.lon=path.lon[y];

			azimuth=(Azimuth(source,temp));

			if (fd!=NULL)
				fprintf(fd,"%.7f, %.7f, %.3f, %.3f, ",path.lat[y], path.lon[y], azimuth, elevation);

			/* If ERP==0, write path loss to alphanumeric
			   output file.  Otherwise, write field strength
			   or received power level (below), as appropriate. */

			if (fd!=NULL && LR.erp==0.0)
				fprintf(fd,"%.2f",loss);

			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			x=(int)rint(10.0*(10.0-elevation));

			if (x>=0 && x<=1000)
			{
				azimuth=rint(azimuth);

				pattern=(double)LR.antenna_pattern[(int)azimuth][x];

				if (pattern!=0.0)
				{
					pattern=20.0*log10(pattern);
					loss-=pattern;
				}
			}

			if (LR.erp!=0.0)
			{
				if (dbm)
				{
					/* dBm is based on EIRP (ERP + 2.14) */

					rxp=LR.erp/(pow(10.0,(loss-2.14)/10.0));

					dBm=10.0*(log10(rxp*1000.0));

					if (fd!=NULL)
						fprintf(fd,"%.3f",dBm);

					/* Scale roughly between 0 and 255 */

					ifs=200+(int)rint(dBm);

					if (ifs<0)
						ifs=0;

					if (ifs>255)
						ifs=255;

					ofs=GetSignal(path.lat[y],path.lon[y]);

					if (ofs>ifs)
						ifs=ofs;

					PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);
				}

				else
				{
					field_strength=(139.4+(20.0*log10(LR.frq_mhz))-loss)+(10.0*log10(LR.erp/1000.0));

					ifs=100+(int)rint(field_strength);

					if (ifs<0)
						ifs=0;

					if (ifs>255)
						ifs=255;

					ofs=GetSignal(path.lat[y],path.lon[y]);

					if (ofs>ifs)
						ifs=ofs;

					PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);

					if (fd!=NULL)
						fprintf(fd,"%.3f",field_strength);
				}
			}

			else
			{
				if (loss>255)
					ifs=255;
				else
					ifs=(int)rint(loss);

				ofs=GetSignal(path.lat[y],path.lon[y]);

				if (ofs<ifs && ofs!=0)
					ifs=ofs;

				PutSignal(path.lat[y],path.lon[y],(unsigned char)ifs);
			}

			if (fd!=NULL)
			{
				if (block)
					fprintf(fd," *");

				fprintf(fd,"\n");
			}

			/* Mark this point as having been analyzed */

			PutMask(path.lat[y],path.lon[y],(GetMask(path.lat[y],path.lon[y])&7)+(mask_value<<3));
		}
	}
}

void PlotLOSMap(struct site source, double altitude)
{
	/* This function performs a 360 degree sweep around the
	   transmitter site (source location), and plots the
	   line-of-sight coverage of the transmitter on the ss
	   generated topographic map based on a receiver located
	   at the specified altitude (in feet AGL).  Results
	   are stored in memory, and written out in the form
	   of a topographic map when the WritePPM() function
	   is later invoked. */

	int y, z, count;
	struct site edge;
	unsigned char symbol[4], x;
	double lat, lon, minwest, maxnorth, th;
	static unsigned char mask_value=1;

	symbol[0]='.';
	symbol[1]='o';
	symbol[2]='O';
	symbol[3]='o';

	count=0;

	//fprintf(stdout,"\nComputing line-of-sight coverage of \"%s\" with an RX antenna\nat %.2f %s AGL",source.name,metric?altitude*METERS_PER_FOOT:altitude,metric?"meters":"feet");

	if (clutter>0.0 && debug)
		fprintf(stdout," and %.2f %s of ground clutter",metric?clutter*METERS_PER_FOOT:clutter,metric?"meters":"feet");

	fprintf(stdout,"...\n\n 0%c to  25%c ",37,37);
	fflush(stdout);

	/* th=pixels/degree divided by 64 loops per
	   progress indicator symbol (.oOo) printed. */

	th=ppd/loops;

	z=(int)(th*ReduceAngle(max_west-min_west));

	minwest=dpp+(double)min_west;
	maxnorth=(double)max_north-dpp;

	for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
	{
		if (lon>=360.0)
			lon-=360.0;

		edge.lat=max_north;
		edge.lon=lon;
		edge.alt=altitude;

		PlotPath(source,edge,mask_value);
		count++;

		if (count==z)
		{
			fprintf(stdout,"%c",symbol[x]);
			fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	count=0;
	fprintf(stdout,"\n25%c to  50%c ",37,37);
	fflush(stdout);

	z=(int)(th*(double)(max_north-min_north));

	for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
	{
		edge.lat=lat;
		edge.lon=min_west;
		edge.alt=altitude;

		PlotPath(source,edge,mask_value);
		count++;

		if (count==z)
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	count=0;
	fprintf(stdout,"\n50%c to  75%c ",37,37);
	fflush(stdout);

	z=(int)(th*ReduceAngle(max_west-min_west));

	for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
	{
		if (lon>=360.0)
			lon-=360.0;

		edge.lat=min_north;
		edge.lon=lon;
		edge.alt=altitude;

		PlotPath(source,edge,mask_value);
		count++;

		if (count==z)
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	count=0;
	fprintf(stdout,"\n75%c to 100%c ",37,37);
	fflush(stdout);

	z=(int)(th*(double)(max_north-min_north));

	for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
	{
		edge.lat=lat;
		edge.lon=max_west;
		edge.alt=altitude;

		PlotPath(source,edge,mask_value);
		count++;

		if (count==z)
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	fprintf(stdout,"\nDone!\n");
	fflush(stdout);

	/* Assign next mask value */

	switch (mask_value)
	{
		case 1:
			mask_value=8;
			break;

		case 8:
			mask_value=16;
			break;

		case 16:
			mask_value=32;
	}
}

void PlotLRMap(struct site source, double altitude, char *plo_filename)
{
	/* This function performs a 360 degree sweep around the
	   transmitter site (source location), and plots the
	   Longley-Rice attenuation on the ss generated
	   topographic map based on a receiver located at
	   the specified altitude (in feet AGL).  Results
	   are stored in memory, and written out in the form
	   of a topographic map when the DoPathLoss() or
	   DoSigStr() functions are later invoked. */

	int y, z, count;
	struct site edge;
	double lat, lon, minwest, maxnorth, th;
	unsigned char x, symbol[4];
	static unsigned char mask_value=1;
	FILE *fd=NULL;

	minwest=dpp+(double)min_west;
	maxnorth=(double)max_north-dpp;

	count=0;
        if (debug){
	fprintf(stdout,"\nComputing Longley-Rice ");
        }

	if (LR.erp==0.0 && debug)
		fprintf(stdout,"path loss");
	else
	{
            if(debug){
		if (dbm)
			fprintf(stdout,"signal power level");
		else
			fprintf(stdout,"field strength");
            }
	}
        if (debug){
	fprintf(stdout," contours of \"%s\"\nout to a radius of %.2f %s with Rx antenna(s) at %.2f %s AGL\n",source.name,metric?max_range*KM_PER_MILE:max_range,metric?"kilometers":"miles",metric?altitude*METERS_PER_FOOT:altitude,metric?"meters":"feet");
        }

	if (clutter>0.0 && debug)
		fprintf(stdout,"\nand %.2f %s of ground clutter",metric?clutter*METERS_PER_FOOT:clutter,metric?"meters":"feet");
        if(debug){
	fprintf(stdout,"...\n\n 0%c to  25%c ",37,37);
	fflush(stdout);
        }
	if (plo_filename[0]!=0)
		fd=fopen(plo_filename,"wb");

	if (fd!=NULL)
	{
		/* Write header information to output file */

		fprintf(fd,"%d, %d\t; max_west, min_west\n%d, %d\t; max_north, min_north\n",max_west, min_west, max_north, min_north);
	}

	/* th=pixels/degree divided by 64 loops per
	   progress indicator symbol (.oOo) printed. */

	th=ppd/loops;

	z=(int)(th*ReduceAngle(max_west-min_west));

	for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
	{
		if (lon>=360.0)
			lon-=360.0;

		edge.lat=max_north;
		edge.lon=lon;
		edge.alt=altitude;

		PlotLRPath(source,edge,mask_value,fd);
		count++;

		if (count==z) 
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	count=0;
        if(debug){
	fprintf(stdout,"\n25%c to  50%c ",37,37);
	fflush(stdout);
        }
	z=(int)(th*(double)(max_north-min_north));

	for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
	{
		edge.lat=lat;
		edge.lon=min_west;
		edge.alt=altitude;

		PlotLRPath(source,edge,mask_value,fd);
		count++;

		if (count==z) 
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	count=0;
        if(debug){
	fprintf(stdout,"\n50%c to  75%c ",37,37);
	fflush(stdout);
        }
	z=(int)(th*ReduceAngle(max_west-min_west));

	for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
	{
		if (lon>=360.0)
			lon-=360.0;

		edge.lat=min_north;
		edge.lon=lon;
		edge.alt=altitude;

		PlotLRPath(source,edge,mask_value,fd);
		count++;
                if (count==z) 
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
		
	}

	count=0;
        if(debug){
	fprintf(stdout,"\n75%c to 100%c ",37,37);
	fflush(stdout);
        }
	z=(int)(th*(double)(max_north-min_north));

	for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
	{
		edge.lat=lat;
		edge.lon=max_west;
		edge.alt=altitude;

		PlotLRPath(source,edge,mask_value,fd);
		count++;

		if (count==z) 
		{
			//fprintf(stdout,"%c",symbol[x]);
			//fflush(stdout);
			count=0;

			if (x==3)
				x=0;
			else
				x++;
		}
	}

	if (fd!=NULL)
		fclose(fd);

	
	if (mask_value<30)
		mask_value++;
}

void LoadSignalColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer=NULL, *s=NULL;
	FILE *fd=NULL;

	for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
		filename[x]=xmtr.filename[x];

	filename[x]='.';
	filename[x+1]='s';
	filename[x+2]='c';
	filename[x+3]='f';
	filename[x+4]=0;

	/* Default values */

	region.level[0]=128;
	region.color[0][0]=255;
	region.color[0][1]=0;
	region.color[0][2]=0;

	region.level[1]=118;
	region.color[1][0]=255;
	region.color[1][1]=165;
	region.color[1][2]=0;

	region.level[2]=108;
	region.color[2][0]=255;
	region.color[2][1]=206;
	region.color[2][2]=0;

	region.level[3]=98;
	region.color[3][0]=255;
	region.color[3][1]=255;
	region.color[3][2]=0;

	region.level[4]=88;
	region.color[4][0]=184;
	region.color[4][1]=255;
	region.color[4][2]=0;

	region.level[5]=78;
	region.color[5][0]=0;
	region.color[5][1]=255;
	region.color[5][2]=0;

	region.level[6]=68;
	region.color[6][0]=0;
	region.color[6][1]=208;
	region.color[6][2]=0;

	region.level[7]=58;
	region.color[7][0]=0;
	region.color[7][1]=196;
	region.color[7][2]=196;

	region.level[8]=48;
	region.color[8][0]=0;
	region.color[8][1]=148;
	region.color[8][2]=255;

	region.level[9]=38;
	region.color[9][0]=80;
	region.color[9][1]=80;
	region.color[9][2]=255;

	region.level[10]=28;
	region.color[10][0]=0;
	region.color[10][1]=38;
	region.color[10][2]=255;

	region.level[11]=18;
	region.color[11][0]=142;
	region.color[11][1]=63;
	region.color[11][2]=255;

	region.level[12]=8;
	region.color[12][0]=140;
	region.color[12][1]=0;
	region.color[12][2]=128;

	region.levels=13;

	fd=fopen(filename,"r");

	if (fd==NULL)
		fd=fopen(filename,"r");

	if (fd==NULL)
	{
		fd=fopen(filename,"w");

		
		for (x=0; x<region.levels; x++)
			fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);
	}

	else
	{
		x=0;
		s=fgets(string,80,fd);

		while (x<128 && feof(fd)==0)
		{
			pointer=strchr(string,';');

			if (pointer!=NULL)
				*pointer=0;

			ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok==4)
			{
				for (y=0; y<4; y++)
				{
					if (val[y]>255)
						val[y]=255;

					if (val[y]<0)
						val[y]=0;
				}

				region.level[x]=val[0];
				region.color[x][0]=val[1];
				region.color[x][1]=val[2];
				region.color[x][2]=val[3];
				x++;
			}

			s=fgets(string,80,fd);
		}

		fclose(fd);
		region.levels=x;
	}
}

void LoadLossColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer=NULL, *s=NULL;
	FILE *fd=NULL;

	for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
		filename[x]=xmtr.filename[x];

	filename[x]='.';
	filename[x+1]='l';
	filename[x+2]='c';
	filename[x+3]='f';
	filename[x+4]=0;

	/* Default values */

	region.level[0]=80;
	region.color[0][0]=255;
	region.color[0][1]=0;
	region.color[0][2]=0;

	region.level[1]=90;
	region.color[1][0]=255;
	region.color[1][1]=128;
	region.color[1][2]=0;

	region.level[2]=100;
	region.color[2][0]=255;
	region.color[2][1]=165;
	region.color[2][2]=0;

	region.level[3]=110;
	region.color[3][0]=255;
	region.color[3][1]=206;
	region.color[3][2]=0;

	region.level[4]=120;
	region.color[4][0]=255;
	region.color[4][1]=255;
	region.color[4][2]=0;

	region.level[5]=130;
	region.color[5][0]=184;
	region.color[5][1]=255;
	region.color[5][2]=0;

	region.level[6]=140;
	region.color[6][0]=0;
	region.color[6][1]=255;
	region.color[6][2]=0;

	region.level[7]=150;
	region.color[7][0]=0;
	region.color[7][1]=208;
	region.color[7][2]=0;

	region.level[8]=160;
	region.color[8][0]=0;
	region.color[8][1]=196;
	region.color[8][2]=196;

	region.level[9]=170;
	region.color[9][0]=0;
	region.color[9][1]=148;
	region.color[9][2]=255;

	region.level[10]=180;
	region.color[10][0]=80;
	region.color[10][1]=80;
	region.color[10][2]=255;

	region.level[11]=190;
	region.color[11][0]=0;
	region.color[11][1]=38;
	region.color[11][2]=255;

	region.level[12]=200;
	region.color[12][0]=142;
	region.color[12][1]=63;
	region.color[12][2]=255;

	region.level[13]=210;
	region.color[13][0]=196;
	region.color[13][1]=54;
	region.color[13][2]=255;

	region.level[14]=220;
	region.color[14][0]=255;
	region.color[14][1]=0;
	region.color[14][2]=255;

	region.level[15]=230;
	region.color[15][0]=255;
	region.color[15][1]=194;
	region.color[15][2]=204;

	region.levels=16;

	fd=fopen(filename,"r");

	if (fd==NULL)
		fd=fopen(filename,"r");

	if (fd==NULL)
	{
		fd=fopen(filename,"w");

		

		for (x=0; x<region.levels; x++)
			fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);
	}

	else
	{
		x=0;
		s=fgets(string,80,fd);

		while (x<128 && feof(fd)==0)
		{
			pointer=strchr(string,';');

			if (pointer!=NULL)
				*pointer=0;

			ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok==4)
			{
				for (y=0; y<4; y++)
				{
					if (val[y]>255)
						val[y]=255;

					if (val[y]<0)
						val[y]=0;
				}

				region.level[x]=val[0];
				region.color[x][0]=val[1];
				region.color[x][1]=val[2];
				region.color[x][2]=val[3];
				x++;
			}

			s=fgets(string,80,fd);
		}

		fclose(fd);
		region.levels=x;
	}
}

void LoadDBMColors(struct site xmtr)
{
	int x, y, ok, val[4];
	char filename[255], string[80], *pointer=NULL, *s=NULL;
	FILE *fd=NULL;

	for (x=0; xmtr.filename[x]!='.' && xmtr.filename[x]!=0 && x<250; x++)
		filename[x]=xmtr.filename[x];

	filename[x]='.';
	filename[x+1]='d';
	filename[x+2]='c';
	filename[x+3]='f';
	filename[x+4]=0;

	/* Default values */

	region.level[0]=0;
	region.color[0][0]=255;
	region.color[0][1]=0;
	region.color[0][2]=0;

	region.level[1]=-10;
	region.color[1][0]=255;
	region.color[1][1]=128;
	region.color[1][2]=0;

	region.level[2]=-20;
	region.color[2][0]=255;
	region.color[2][1]=165;
	region.color[2][2]=0;

	region.level[3]=-30;
	region.color[3][0]=255;
	region.color[3][1]=206;
	region.color[3][2]=0;

	region.level[4]=-40;
	region.color[4][0]=255;
	region.color[4][1]=255;
	region.color[4][2]=0;

	region.level[5]=-50;
	region.color[5][0]=184;
	region.color[5][1]=255;
	region.color[5][2]=0;

	region.level[6]=-60;
	region.color[6][0]=0;
	region.color[6][1]=255;
	region.color[6][2]=0;

	region.level[7]=-70;
	region.color[7][0]=0;
	region.color[7][1]=208;
	region.color[7][2]=0;

	region.level[8]=-80;
	region.color[8][0]=0;
	region.color[8][1]=196;
	region.color[8][2]=196;

	region.level[9]=-90;
	region.color[9][0]=0;
	region.color[9][1]=148;
	region.color[9][2]=255;

	region.level[10]=-100;
	region.color[10][0]=80;
	region.color[10][1]=80;
	region.color[10][2]=255;

	region.level[11]=-110;
	region.color[11][0]=0;
	region.color[11][1]=38;
	region.color[11][2]=255;

	region.level[12]=-120;
	region.color[12][0]=142;
	region.color[12][1]=63;
	region.color[12][2]=255;

	region.level[13]=-130;
	region.color[13][0]=196;
	region.color[13][1]=54;
	region.color[13][2]=255;

	region.level[14]=-140;
	region.color[14][0]=255;
	region.color[14][1]=0;
	region.color[14][2]=255;

	region.level[15]=-150;
	region.color[15][0]=255;
	region.color[15][1]=194;
	region.color[15][2]=204;

	region.levels=16;

	fd=fopen(filename,"r");

	if (fd==NULL)
		fd=fopen(filename,"r");

	if (fd==NULL)
	{
		fd=fopen(filename,"w");

		
		for (x=0; x<region.levels; x++)
			fprintf(fd,"%+4d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

		fclose(fd);
	}

	else
	{
		x=0;
		s=fgets(string,80,fd);

		while (x<128 && feof(fd)==0)
		{
			pointer=strchr(string,';');

			if (pointer!=NULL)
				*pointer=0;

			ok=sscanf(string,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

			if (ok==4)
			{
				if (val[0]<-200)
					val[0]=-200;

				if (val[0]>+40)
					val[0]=+40;

				region.level[x]=val[0];

				for (y=1; y<4; y++)
				{
					if (val[y]>255)
						val[y]=255;

					if (val[y]<0)
						val[y]=0;
				}

				region.color[x][0]=val[1];
				region.color[x][1]=val[2];
				region.color[x][2]=val[3];
				x++;
			}

			s=fgets(string,80,fd);
		}

		fclose(fd);
		region.levels=x;
	}
}


void DoPathLoss(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the content of flags held in the mask[][]
	   array (only).  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	char mapfile[255], geofile[255], kmlfile[255];
	unsigned width, height, red, green, blue, terrain=0;
	unsigned char found, mask, cityorcounty;
	int indx, x, y, z, x0, y0, loss, match;
	double lat, lon, conversion, one_over_gamma,minwest;
	FILE *fd;

	one_over_gamma=1.0/GAMMA;
	conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

	width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
	height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

	LoadLossColors(xmtr[0]);

	if (filename[0]==0)
	{
		strncpy(filename, xmtr[0].filename,254);
		filename[strlen(filename)-4]=0;  /* Remove .qth */
	}

	y=strlen(filename);

	if (y>4)
	{
		if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
			y-=4;
	}

	for (x=0; x<y; x++)
	{
		mapfile[x]=filename[x];
		geofile[x]=filename[x];
		kmlfile[x]=filename[x];
	}

	mapfile[x]='.';
	geofile[x]='.';
	kmlfile[x]='.';
	mapfile[x+1]='p';
	geofile[x+1]='g';
	kmlfile[x+1]='k';
	mapfile[x+2]='p';
	geofile[x+2]='e';
	kmlfile[x+2]='m';
	mapfile[x+3]='m';
	geofile[x+3]='o';
	kmlfile[x+3]='l';
	mapfile[x+4]=0;
	geofile[x+4]=0;
	kmlfile[x+4]=0;

	minwest=((double)min_west)+dpp;

	if (minwest>360.0)
		minwest-=360.0;

	north=(double)max_north-dpp;

	if (kml || geo)
		south=(double)min_north;	/* No bottom legend */
	else
		south=(double)min_north-(30.0/ppd); /* 30 pixels for bottom legend */

	east=(minwest<180.0?-minwest:360.0-min_west);
	west=(double)(max_west<180?-max_west:360-max_west);

	
	// WriteKML()
        //writeKML(xmtr,filename);

	fd=fopen(mapfile,"wb");

	fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
	fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
	fflush(stdout);
        }
	for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
	{
		for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
		{
			if (lon<0.0)
				lon+=360.0;

			for (indx=0, found=0; indx<MAXPAGES && found==0;)
			{
				x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
				y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

				if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
					found=1;
				else
					indx++;
			}
                        
                        
                        
			if (found)
			{
				mask=dem[indx].mask[x0][y0];
				loss=(dem[indx].signal[x0][y0]);
				cityorcounty=0;

				match=255;

				red=0;
				green=0;
				blue=0;

				if (loss<=region.level[0])
					match=0;
				else
				{
					for (z=1; (z<region.levels && match==255); z++)
					{
						if (loss>=region.level[z-1] && loss<region.level[z])
							match=z;
					}
				}

				if (match<region.levels)
				{
					red=region.color[match][0];
					green=region.color[match][1];
					blue=region.color[match][2];
				}

	 			if (mask&2)
				{
					/* Text Labels: Red or otherwise */

					if (red>=180 && green<=75 && blue<=75 && loss==0)
                                                fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
                                        else
                                                fprintf(fd,"%c%c%c",255,0,0);

                                        cityorcounty=1;
				}

				else if (mask&4)
				{
					/* County Boundaries: Black */

					fprintf(fd,"%c%c%c",0,0,0);

					cityorcounty=1;
				}

				if (cityorcounty==0)
				{
					if (loss==0 || (contour_threshold!=0 && loss>abs(contour_threshold)))
					{
						if (ngs)  /* No terrain */
							fprintf(fd,"%c%c%c",255,255,255);
						else
						{
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0]==0)
								fprintf(fd,"%c%c%c",0,0,170);
							else
							{
								terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
								fprintf(fd,"%c%c%c",terrain,terrain,terrain);
							}
						}
					}

					else
					{
						/* Plot path loss in color */

						if (red!=0 || green!=0 || blue!=0)
							fprintf(fd,"%c%c%c",red,green,blue);

						else  /* terrain / sea-level */
						{
							if (dem[indx].data[x0][y0]==0)
								fprintf(fd,"%c%c%c",0,0,170);
							else
							{
								/* Elevation: Greyscale */
								terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
								fprintf(fd,"%c%c%c",terrain,terrain,terrain);
							}
						}
					}
				}
			}

			else
			{
				/* We should never get here, but if */
				/* we do, display the region as black */

				fprintf(fd,"%c%c%c",0,0,0);
			}
		}
	}



	fclose(fd);
        if(debug){
	fprintf(stdout,"Done!\n");
	fflush(stdout);
        }
}

void DoSigStr(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the signal strength values held in the
	   signal[][] array.  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	char mapfile[255], geofile[255], kmlfile[255];
	unsigned width, height, terrain, red, green, blue;
	unsigned char found, mask, cityorcounty;
	int indx, x, y, z=1, x0, y0, signal,match;
	double conversion, one_over_gamma, lat, lon, minwest;
	FILE *fd;

	one_over_gamma=1.0/GAMMA;
	conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

	width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
	height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

	LoadSignalColors(xmtr[0]);

	if (filename[0]==0)
	{
		strncpy(filename, xmtr[0].filename,254);
		filename[strlen(filename)-4]=0;  /* Remove .qth */
	}

	y=strlen(filename);

	if (y>4)
	{
		if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
			y-=4;
	}

	for (x=0; x<y; x++)
	{
		mapfile[x]=filename[x];
		geofile[x]=filename[x];
		kmlfile[x]=filename[x];
	}

	mapfile[x]='.';
	geofile[x]='.';
	kmlfile[x]='.';
	mapfile[x+1]='p';
	geofile[x+1]='g';
	kmlfile[x+1]='k';
	mapfile[x+2]='p';
	geofile[x+2]='e';
	kmlfile[x+2]='m';
	mapfile[x+3]='m';
	geofile[x+3]='o';
	kmlfile[x+3]='l';
	mapfile[x+4]=0;
	geofile[x+4]=0;
	kmlfile[x+4]=0;

	minwest=((double)min_west)+dpp;

	if (minwest>360.0)
		minwest-=360.0;

	north=(double)max_north-dpp;

	
	south=(double)min_north;	/* No bottom legend */
	
	east=(minwest<180.0?-minwest:360.0-min_west);
	west=(double)(max_west<180?-max_west:360-max_west);

	// WriteKML()
        //writeKML(xmtr,filename);

	fd=fopen(mapfile,"wb");

	fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
	fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
	fflush(stdout);
        }
	for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
	{
		for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
		{
			if (lon<0.0)
				lon+=360.0;

			for (indx=0, found=0; indx<MAXPAGES && found==0;)
			{
				x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
				y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

				if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
					found=1;
				else
					indx++;
			}

			if (found)
			{
				mask=dem[indx].mask[x0][y0];
				signal=(dem[indx].signal[x0][y0])-100;
				cityorcounty=0;

				match=255;

				red=0;
				green=0;
				blue=0;

				if (signal>=region.level[0])
					match=0;
				else
				{
					for (z=1; (z<region.levels && match==255); z++)
					{
						if (signal<region.level[z-1] && signal>=region.level[z])
							match=z;
					}
				}

				if (match<region.levels)
				{
					red=region.color[match][0];
					green=region.color[match][1];
					blue=region.color[match][2];
				}

	 			if (mask&2)
				{
					/* Text Labels: Red or otherwise */

					if (red>=180 && green<=75 && blue<=75)
						fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
					else
						fprintf(fd,"%c%c%c",255,0,0);

					cityorcounty=1;
				}

				else if (mask&4)
				{
					/* County Boundaries: Black */

					fprintf(fd,"%c%c%c",0,0,0);

					cityorcounty=1;
				}

				if (cityorcounty==0)
				{
					if (contour_threshold!=0 && signal<contour_threshold)
					{
						if (ngs)
							fprintf(fd,"%c%c%c",255,255,255);
						else
						{
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0]==0)
								fprintf(fd,"%c%c%c",0,0,170);
							else
							{
								terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
								fprintf(fd,"%c%c%c",terrain,terrain,terrain);
							}
						}
					}

					else
					{
						/* Plot field strength regions in color */

						if (red!=0 || green!=0 || blue!=0)
							fprintf(fd,"%c%c%c",red,green,blue);

						else  /* terrain / sea-level */
						{
							if (ngs)
								fprintf(fd,"%c%c%c",255,255,255);
							else
							{
								if (dem[indx].data[x0][y0]==0)
									fprintf(fd,"%c%c%c",0,0,170);
								else
								{
									/* Elevation: Greyscale */
									terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
									fprintf(fd,"%c%c%c",terrain,terrain,terrain);
								}
							}
						}
					}
				}
			}

			else
			{
				/* We should never get here, but if */
				/* we do, display the region as black */

				fprintf(fd,"%c%c%c",0,0,0);
			}
		}
	}


	fclose(fd);
        if(debug){
	fprintf(stdout,"Done!\n");
	fflush(stdout);
        }
}

void DoRxdPwr(char *filename, unsigned char geo, unsigned char kml, unsigned char ngs, struct site *xmtr, unsigned char txsites)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the signal power level values held in the
	   signal[][] array.  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	char mapfile[255], geofile[255], kmlfile[255];
	unsigned width, height, terrain, red, green, blue;
	unsigned char found, mask, cityorcounty;
	int indx, x, y, z=1, x0, y0, dBm, match;
	double conversion, one_over_gamma, lat, lon, minwest;
	FILE *fd;

	one_over_gamma=1.0/GAMMA;
	conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

	width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
	height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

	LoadDBMColors(xmtr[0]);

	if (filename[0]==0)
	{
		strncpy(filename, xmtr[0].filename,254);
		filename[strlen(filename)-4]=0;  /* Remove .qth */
	}

	y=strlen(filename);

	if (y>4)
	{
		if (filename[y-1]=='m' && filename[y-2]=='p' && filename[y-3]=='p' && filename[y-4]=='.')
			y-=4;
	}

	for (x=0; x<y; x++)
	{
		mapfile[x]=filename[x];
		geofile[x]=filename[x];
		kmlfile[x]=filename[x];
	}

	mapfile[x]='.';
	geofile[x]='.';
	kmlfile[x]='.';
	mapfile[x+1]='p';
	geofile[x+1]='g';
	kmlfile[x+1]='k';
	mapfile[x+2]='p';
	geofile[x+2]='e';
	kmlfile[x+2]='m';
	mapfile[x+3]='m';
	geofile[x+3]='o';
	kmlfile[x+3]='l';
	mapfile[x+4]=0;
	geofile[x+4]=0;
	kmlfile[x+4]=0;

	minwest=((double)min_west)+dpp;

	if (minwest>360.0)
		minwest-=360.0;

	north=(double)max_north-dpp;

	
	south=(double)min_north;	/* No bottom legend */
	

	east=(minwest<180.0?-minwest:360.0-min_west);
	west=(double)(max_west<180?-max_west:360-max_west);

	fd=fopen(mapfile,"wb");

	fprintf(fd,"P6\n%u %u\n255\n",width,(kml?height:height+30));
        if(debug){
	fprintf(stdout,"\nWriting \"%s\" (%ux%u pixmap image)... ",mapfile,width,(kml?height:height+30));
	fflush(stdout);
        }
        // WriteKML()
        //writeKML(xmtr,filename);
        
	for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
	{
		for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
		{
			if (lon<0.0)
				lon+=360.0;

			for (indx=0, found=0; indx<MAXPAGES && found==0;)
			{
				x0=(int)rint(ppd*(lat-(double)dem[indx].min_north));
				y0=mpi-(int)rint(ppd*(LonDiff((double)dem[indx].max_west,lon)));

				if (x0>=0 && x0<=mpi && y0>=0 && y0<=mpi)
					found=1;
				else
					indx++;
			}

			if (found)
			{
				mask=dem[indx].mask[x0][y0];
				dBm=(dem[indx].signal[x0][y0])-200;
				cityorcounty=0;

				match=255;

				red=0;
				green=0;
				blue=0;

				if (dBm>=region.level[0])
					match=0;
				else
				{
					for (z=1; (z<region.levels && match==255); z++)
					{
						if (dBm<region.level[z-1] && dBm>=region.level[z])
							match=z;
					}
				}

				if (match<region.levels)
				{
					red=region.color[match][0];
					green=region.color[match][1];
					blue=region.color[match][2];
				}

	 			if (mask&2)
				{
					/* Text Labels: Red or otherwise */

					if (red>=180 && green<=75 && blue<=75 && dBm!=0)
						fprintf(fd,"%c%c%c",255^red,255^green,255^blue);
					else
						fprintf(fd,"%c%c%c",255,0,0);

					cityorcounty=1;
				}

				else if (mask&4)
				{
					/* County Boundaries: Black */

					fprintf(fd,"%c%c%c",0,0,0);

					cityorcounty=1;
				}

				if (cityorcounty==0)
				{
					if (contour_threshold!=0 && dBm<contour_threshold)
					{
						if (ngs) /* No terrain */
							fprintf(fd,"%c%c%c",255,255,255);
						else
						{
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0]==0)
								fprintf(fd,"%c%c%c",0,0,170);
							else
							{
								terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
								fprintf(fd,"%c%c%c",terrain,terrain,terrain);
							}
						}
					}

					else
					{
						/* Plot signal power level regions in color */

						if (red!=0 || green!=0 || blue!=0)
							fprintf(fd,"%c%c%c",red,green,blue);

						else  /* terrain / sea-level */
						{
							if (ngs)
								fprintf(fd,"%c%c%c",255,255,255);
							else
							{
								if (dem[indx].data[x0][y0]==0)
									fprintf(fd,"%c%c%c",0,0,170);
								else
								{
									/* Elevation: Greyscale */
									terrain=(unsigned)(0.5+pow((double)(dem[indx].data[x0][y0]-min_elevation),one_over_gamma)*conversion);
									fprintf(fd,"%c%c%c",terrain,terrain,terrain);
								}
							}
						}
					}
				}
			}

			else
			{
				/* We should never get here, but if */
				/* we do, display the region as black */

				fprintf(fd,"%c%c%c",0,0,0);
			}
		}
	}



	fclose(fd);
        if(debug){
	fprintf(stdout,"Done!\n");
	fflush(stdout);
        }
}

void LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat)
{
	/* This function loads the SDF files required
	   to cover the limits of the region specified. */

	int x, y, width, ymin, ymax;

	width=ReduceAngle(max_lon-min_lon);

	if ((max_lon-min_lon)<=180.0)
	{
		for (y=0; y<=width; y++)
			for (x=min_lat; x<=max_lat; x++)
			{
				ymin=(int)(min_lon+(double)y);

				while (ymin<0)
					ymin+=360;

				while (ymin>=360)
					ymin-=360;

				ymax=ymin+1;

				while (ymax<0)
					ymax+=360;

				while (ymax>=360)
					ymax-=360;

				if (ippd==3600)
					snprintf(string,19,"%d:%d:%d:%d-hd",x, x+1, ymin, ymax);
				else
					snprintf(string,16,"%d:%d:%d:%d",x, x+1, ymin, ymax);
                                

				LoadSDF(string);
			}
	}

	else
	{
		for (y=0; y<=width; y++)
			for (x=min_lat; x<=max_lat; x++)
			{
				ymin=max_lon+y;

				while (ymin<0)
					ymin+=360;

				while (ymin>=360)
					ymin-=360;

				ymax=ymin+1;

				while (ymax<0)
					ymax+=360;

				while (ymax>=360)
					ymax-=360;

				if (ippd==3600)
					snprintf(string,19,"%d:%d:%d:%d-hd",x, x+1, ymin, ymax);
				else
					snprintf(string,16,"%d:%d:%d:%d",x, x+1, ymin, ymax);
				LoadSDF(string);
			}
	}
}



void LoadUDT(char *filename)
{
	/* This function reads a file containing User-Defined Terrain
	   features for their addition to the digital elevation model
	   data used by SPLAT!.  Elevations in the UDT file are evaluated
	   and then copied into a temporary file under /tmp.  Then the
	   contents of the temp file are scanned, and if found to be unique,
	   are added to the ground elevations described by the digital
	   elevation data already loaded into memory. */

	int	i, x, y, z, ypix, xpix, tempxpix, tempypix, fd=0, n=0, pixelfound=0;
	char	input[80], str[3][80], tempname[15], *pointer=NULL, *s=NULL;
	double	latitude, longitude, height, tempheight;
	FILE	*fd1=NULL, *fd2=NULL;

	strcpy(tempname,"/tmp/XXXXXX\0");

	fd1=fopen(filename,"r");

	if (fd1!=NULL)
	{
		fd=mkstemp(tempname);
		fd2=fopen(tempname,"w");

		s=fgets(input,78,fd1);

		pointer=strchr(input,';');

		if (pointer!=NULL)
			*pointer=0;

		
		while (feof(fd1)==0)
		{
			/* Parse line for latitude, longitude, height */

			for (x=0, y=0, z=0; x<78 && input[x]!=0 && z<3; x++)
			{
				if (input[x]!=',' && y<78)
				{
					str[z][y]=input[x];
					y++;
				}

				else
				{
					str[z][y]=0;
					z++;
					y=0;
				}
			}

			latitude=ReadBearing(str[0]);
			longitude=ReadBearing(str[1]);

			if (longitude<0.0)
				longitude+=360; 
			
			/* Remove <CR> and/or <LF> from antenna height string */

			for (i=0; str[2][i]!=13 && str[2][i]!=10 && str[2][i]!=0; i++);

			str[2][i]=0;

			/* The terrain feature may be expressed in either
			   feet or meters.  If the letter 'M' or 'm' is
			   discovered in the string, then this is an
			   indication that the value given is expressed
			   in meters.  Otherwise the height is interpreted
			   as being expressed in feet.  */

			for (i=0; str[2][i]!='M' && str[2][i]!='m' && str[2][i]!=0 && i<48; i++);

			if (str[2][i]=='M' || str[2][i]=='m')
			{
				str[2][i]=0;
				height=rint(atof(str[2]));
			}

			else
			{
				str[2][i]=0;
				height=rint(METERS_PER_FOOT*atof(str[2]));
			}

			if (height>0.0)
				fprintf(fd2,"%d, %d, %f\n",(int)rint(latitude/dpp), (int)rint(longitude/dpp), height);
        

			s=fgets(input,78,fd1);

			pointer=strchr(input,';');

			if (pointer!=NULL)
				*pointer=0;
		}

		fclose(fd1);
		fclose(fd2);
		close(fd);

		
		fd1=fopen(tempname,"r");
		fd2=fopen(tempname,"r");

		y=0;

		n=fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);

		do
		{
			x=0;
			z=0;

			n=fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);

			do
			{
				if (x>y && xpix==tempxpix && ypix==tempypix)
				{
						z=1;  /* Dupe! */

						if (tempheight>height)
							height=tempheight;
				}

				else
				{
					n=fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);
					x++;
				}

			} while (feof(fd2)==0 && z==0);

			if (z==0)  /* No duplicate found */
		
				//fprintf(stdout,"%lf, %lf \n",xpix*dpp, ypix*dpp);
				fflush(stdout);
				pixelfound = AddElevation(xpix*dpp, ypix*dpp, height);
				//fprintf(stdout,"%d \n",pixelfound);
				fflush(stdout);
				
			n=fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);
			y++;

			rewind(fd2);

		} while (feof(fd1)==0);

		fclose(fd1);
		fclose(fd2);
		unlink(tempname);
	}

	//else
		//fprintf(stderr,"\n*** ERROR: \"%s\": not found!",filename);

	//fprintf(stdout,"\n");
}


int main(int argc, char *argv[])
{
	int		x, y, z=0, min_lat, min_lon, max_lat, max_lon,
			rxlat, rxlon, txlat, txlon, west_min, west_max,
			north_min, north_max;

	unsigned char	LRmap=0, map=0,txsites=0,
			topomap=0, geo=0, kml=0, area_mode=0, max_txsites, ngs=0;

	char		mapfile[255], elevation_file[255], longley_file[255], terrain_file[255],
			string[255], rxfile[255],txfile[255], udt_file[255], rxsite=0, ani_filename[255],
			ano_filename[255];

	double		altitude=0.0, altitudeLR=0.0, tx_range=0.0,
			rx_range=0.0, deg_range=0.0, deg_limit=0.0,
			deg_range_lon;

	struct		site tx_site[32], rx_site;

	

	strncpy(ss_version,"1.3.4\0",6);
	strncpy(ss_name,"Signal Server\0",14);
	
	if (argc==1)
	{
		fprintf(stdout,"\n\t\t -- %s %s options --\n\n",ss_name, ss_version);
                fprintf(stdout,"       -d Directory containing .sdf tiles\n");
                fprintf(stdout,"     -lat Tx Latitude (decimal degrees)\n");
                fprintf(stdout,"     -lon Tx Longitude (decimal degrees) Positive 0-360 \n");
                fprintf(stdout,"     -txh Tx Height (above ground)\n");
                fprintf(stdout,"       -f Tx Frequency (MHz)\n");
                fprintf(stdout,"     -erp Tx Effective Radiated Power (Watts)\n");
		            fprintf(stdout,"     -rxh Rx Height(s) (optional. Default=0.1)\n");
                fprintf(stdout,"      -rt Rx Threshold (dB / dBm / dBuV/m)\n");
                fprintf(stdout,"      -hp Horizontal Polarisation (default=vertical)\n");
		            fprintf(stdout,"      -gc Ground clutter (feet/meters)\n");
                fprintf(stdout,"      -udt User defined terrain filename\n");
	             	fprintf(stdout,"     -dbm Plot Rxd signal power instead of field strength\n");
	             	fprintf(stdout,"       -m Metric units of measurement\n");
                fprintf(stdout,"      -te Terrain code 1-6 (optional)\n");
                fprintf(stdout,"      -terdic Terrain dielectric value 2-80 (optional)\n");
	             	fprintf(stdout,"      -tercon Terrain conductivity 0.01-0.0001 (optional)\n");
                fprintf(stdout,"      -cl Climate code 1-6 (optional)\n");
                fprintf(stdout,"       -o Filename. Required. \n");
                fprintf(stdout,"       -R Radius (miles/kilometers)\n");
                fprintf(stdout,"     -res Pixels per degree. 300/600/1200(default)/3600 (optional)\n");
                fprintf(stdout,"     -t Terrain background\n");
		fprintf(stdout,"     -dbg Debug\n\n");

 

		fflush(stdout);

		return 1;
	}

	y=argc-1;

	kml=0;
	geo=0;
	dbm=0;
	gpsav=0;
	metric=0;
	rxfile[0]=0;
	txfile[0]=0;
	string[0]=0;
	mapfile[0]=0;
	clutter=0.0;
	forced_erp=-1.0;
	forced_freq=0.0;
	elevation_file[0]=0;
	terrain_file[0]=0;
	sdf_path[0]=0;
	udt_file[0]=0;
	path.length=0;
	max_txsites=30;
	fzone_clearance=0.6;
	contour_threshold=0;
	rx_site.lat=91.0;
	rx_site.lon=361.0;
	longley_file[0]=0;
	ano_filename[0]=0;
	ani_filename[0]=0;
	earthradius=EARTHRADIUS;
        max_range=1.0;
                

        lat=0;
        lon=0;
        txh=0;
        ngs=1; // no terrain background
        kml=1;
        
        map=1;
	LRmap=1;
	area_mode=1;
        ippd=IPPD; // default resolution

        sscanf("0.1","%lf",&altitudeLR);

                // Defaults
                LR.eps_dielect=15.0; // Farmland
		LR.sgm_conductivity=0.005; // Farmland
		LR.eno_ns_surfref=301.0; 
		LR.frq_mhz=19.0; // Deliberately too low
		LR.radio_climate=5; // continental
		LR.pol=1; // vert
		LR.conf=0.50; 
		LR.rel=0.50; 
		LR.erp=0.0; // will default to Path Loss

		tx_site[0].lat=91.0;
		tx_site[0].lon=361.0;
	

	for (x=0; x<MAXPAGES; x++)
	{
		dem[x].min_el=32768;
		dem[x].max_el=-32768;
		dem[x].min_north=90;
		dem[x].max_north=-90;
		dem[x].min_west=360;
		dem[x].max_west=-1;
	}

	/* Scan for command line arguments */

	for (x=1; x<=y; x++)
	{

            if (strcmp(argv[x],"-res")==0)
		{
                z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
                                           
                             sscanf(argv[z],"%d",&ippd);
                             
                            switch (ippd)
                                {
                                case 300: 
                                    MAXRAD=300;
                                    jgets=3;
                                break;

                                case 600: 
                                    MAXRAD=150;
                                    jgets=1;
                                break;
				
				case 3600: 
                                    MAXRAD=100;
				    ippd=3600;
                                    jgets=0;
                                break;

                                
                                default:
                                    MAXRAD=100;
                                    ippd=1200;
                                    jgets=0;
                                    break;
                        }
                        }
            }
		if (strcmp(argv[x],"-R")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&max_range);
				
			}
		}


		if (strcmp(argv[x],"-gc")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&clutter);

				if (clutter<0.0)
					clutter=0.0;
			}
		}


    
    
    
		if (strcmp(argv[x],"-o")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
                        {
				 strncpy(mapfile,argv[z],253);
                                 strncpy (tx_site[0].name,"Tx",2);
                                 strncpy (tx_site[0].filename,argv[z],253);
                                 LoadPAT(argv[z]);
                                                              
                        }
			map=1;
		}

            
		if (strcmp(argv[x],"-rt")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0]) /* A minus argument is legal here */
				sscanf(argv[z],"%d",&contour_threshold);
		}

		if (strcmp(argv[x],"-m")==0)
                {
			metric=1;
                        
                }

		if (strcmp(argv[x],"-t")==0)
		{
		ngs=0; // greyscale background
		}

		if (strcmp(argv[x],"-dbm")==0)
			dbm=1;

		if (strcmp(argv[x],"-d")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(sdf_path,argv[z],253);
		}

		if (strcmp(argv[x],"-lat")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0]) // Now allowing for southern hemisphere!
			{
					tx_site[0].lat = ReadBearing(argv[z]);
			}



		}
                if (strcmp(argv[x],"-lon")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0]) /* A minus argument is legal here */
			{
			
				tx_site[0].lon = ReadBearing(argv[z]);
				
				tx_site[0].lon*=-1;

				if (tx_site[0].lon<0.0)
					tx_site[0].lon+=360.0;
                        }



		}
                 if (strcmp(argv[x],"-txh")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
                            sscanf(argv[z],"%f",&tx_site[0].alt);
        		}
                    txsites=1;
		}


		if (strcmp(argv[x],"-rxh")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&altitudeLR);

			}
		}


		if (strcmp(argv[x],"-f")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&LR.frq_mhz);
			}
		}

		if (strcmp(argv[x],"-erp")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				sscanf(argv[z],"%lf",&LR.erp);
			}
		}
               

                if (strcmp(argv[x],"-cl")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{

                            sscanf(argv[z],"%d",&LR.radio_climate);

			}
		}
                 if (strcmp(argv[x],"-te")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{

                            sscanf(argv[z],"%d",&ter);

                          switch (ter)
                                {
                                    case 1: // Water
                                    terdic = 80;
                                    tercon = 0.010;
                                    break;

                                    case 2: // Marsh
                                    terdic = 12;
                                    tercon = 0.007;
                                    break;

                               	    case 3: // Farmland
                                    terdic = 15;
                                    tercon = 0.005;
                                    break;

                                    case 4: //Mountain
                                    terdic = 13;
                                    tercon = 0.002;
                                    break;
                                     case 5: //Desert
                                    terdic=13;
                                    tercon=0.002;
                                    break;
                                    case 6: //Urban
                                    terdic=5;
                                    tercon=0.001;
                                    break;
                        }

                          LR.eps_dielect = terdic;
                          LR.sgm_conductivity = tercon;


			}
		}

		if (strcmp(argv[x],"-terdic")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{

                           sscanf(argv[z],"%lf",&terdic);

                          LR.eps_dielect = terdic;
                         
			}
		}

		if (strcmp(argv[x],"-tercon")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0] && argv[z][0]!='-')
			{

                           sscanf(argv[z],"%lf",&tercon);

                          LR.sgm_conductivity = tercon;
                         
			}
		}

		if (strcmp(argv[x],"-hp")==0)
		{
			// Horizontal polarisation (0)
                    LR.pol = 0;
		}

                if (strcmp(argv[x],"-dbg")==0)
		{
                   debug=1;
		}

    
                ppd=(double)ippd;	/* pixels per degree (double)  */
                dpp=1.0/ppd;		/* degrees per pixel */
                mpi=ippd-1;		/* maximum pixel index per degree */

                
    
    /*UDT*/
   
    if (strcmp(argv[x],"-udt")==0)
		{
			z=x+1;

			if (z<=y && argv[z][0])
			{
           strncpy(udt_file,argv[z],253);
			}
		}
    

	}

        /* ERROR DETECTION */
                 if (tx_site[0].lat > 90 || tx_site[0].lat < -90)
		{
                      fprintf(stdout,"ERROR: Either the lat was missing or out of range!");
                      exit(0);

                 }
                  if (tx_site[0].lon > 360 || tx_site[0].lon < 0)
		{
                      fprintf(stdout,"ERROR: Either the lon was missing or out of range!");
                      exit(0);

                 }

                if (LR.frq_mhz < 20 || LR.frq_mhz > 20000)
                {
                    fprintf(stdout,"ERROR: Either the Frequency was missing or out of range!");
                      exit(0);
                }
                if (LR.erp>2000000)
                {
                     fprintf(stdout,"ERROR: Power was out of range!");
                      exit(0);

                }
		if (LR.eps_dielect > 80 || LR.eps_dielect < 0.1)
                {
                     fprintf(stdout,"ERROR: Ground Dielectric value out of range!");
                      exit(0);

                }
		if (LR.sgm_conductivity > 0.01 || LR.sgm_conductivity < 0.000001)
                {
                     fprintf(stdout,"ERROR: Ground conductivity out of range!");
                      exit(0);

                }

                if (tx_site[0].alt < 0 || tx_site[0].alt > 60000)
                {
                     fprintf(stdout,"ERROR: Tx altitude above ground was too high!");
                      exit(0);
                }
                if (altitudeLR < 0 || altitudeLR > 60000)
                {
                     fprintf(stdout,"ERROR: Rx altitude above ground was too high!");
                      exit(0);
                }


                if (ippd < 300 || ippd > 3600){
                    fprintf(stdout,"ERROR: resolution out of range!");
                    exit(0);
                }

                if(contour_threshold < -200 || contour_threshold > 200)
                {
                    fprintf(stdout,"ERROR: Receiver threshold out of range (-200 / +200)");
                    exit(0);
                }


	 /* ERROR DETECTION COMPLETE */

             
                
	if (metric)
	{
		altitudeLR/=METERS_PER_FOOT;	/* RXH meters --> feet */
		max_range/=KM_PER_MILE;		/* RAD kilometers --> miles */
		//altitude/=METERS_PER_FOOT;	
		tx_site[0].alt/=METERS_PER_FOOT;	/* TXH meters --> feet */
		clutter/=METERS_PER_FOOT;		/* CLH meters --> feet */
		
	}

	/* Ensure a trailing '/' is present in sdf_path */

	if (sdf_path[0])
	{
		x=strlen(sdf_path);

		if (sdf_path[x-1]!='/' && x!=0)
		{
			sdf_path[x]='/';
			sdf_path[x+1]=0;
		}
	}

	x=0;
	y=0;

	min_lat=70;
	max_lat=-70;

	min_lon=(int)floor(tx_site[0].lon);
	max_lon=(int)floor(tx_site[0].lon);

	
		txlat=(int)floor(tx_site[0].lat);
		txlon=(int)floor(tx_site[0].lon);

		if (txlat<min_lat)
			min_lat=txlat;

		if (txlat>max_lat)
			max_lat=txlat;

		if (LonDiff(txlon,min_lon)<0.0)
			min_lon=txlon;

		if (LonDiff(txlon,max_lon)>=0.0)
			max_lon=txlon;
	

	if (rxsite)
	{
		rxlat=(int)floor(rx_site.lat);
		rxlon=(int)floor(rx_site.lon);

		if (rxlat<min_lat)
			min_lat=rxlat;

		if (rxlat>max_lat)
			max_lat=rxlat;

		if (LonDiff(rxlon,min_lon)<0.0)
			min_lon=rxlon;

		if (LonDiff(rxlon,max_lon)>=0.0)
			max_lon=rxlon;
	}

	/* Load the required SDF files */

	LoadTopoData(max_lon, min_lon, max_lat, min_lat);

	if (area_mode || topomap)
	{
		for (z=0; z<txsites && z<max_txsites; z++)
		{
			/* "Ball park" estimates used to load any additional
			   SDF files required to conduct this analysis. */

			tx_range=sqrt(1.5*(tx_site[z].alt+GetElevation(tx_site[z])));

			if (LRmap)
				rx_range=sqrt(1.5*altitudeLR);
			else
				rx_range=sqrt(1.5*altitude);

			/* deg_range determines the maximum
			   amount of topo data we read */

			deg_range=(tx_range+rx_range)/57.0;

			/* max_range regulates the size of the
			   analysis.  A small, non-zero amount can
			   be used to shrink the size of the analysis
			   and limit the amount of topo data read by
			   ss  A large number will increase the
			   width of the analysis and the size of
			   the map. */

			if (max_range==0.0)
				max_range=tx_range+rx_range;

			deg_range=max_range/57.0;

			// No more than 8 degs
			deg_limit=3.5;


			if (fabs(tx_site[z].lat)<70.0)
				deg_range_lon=deg_range/cos(DEG2RAD*tx_site[z].lat);
			else
				deg_range_lon=deg_range/cos(DEG2RAD*70.0);

			/* Correct for squares in degrees not being square in miles */

			if (deg_range>deg_limit)
				deg_range=deg_limit;

			if (deg_range_lon>deg_limit)
				deg_range_lon=deg_limit;

			north_min=(int)floor(tx_site[z].lat-deg_range);
			north_max=(int)floor(tx_site[z].lat+deg_range);

			west_min=(int)floor(tx_site[z].lon-deg_range_lon);

			while (west_min<0)
				west_min+=360;

			while (west_min>=360)
				west_min-=360;

			west_max=(int)floor(tx_site[z].lon+deg_range_lon);

			while (west_max<0)
				west_max+=360;

			while (west_max>=360)
				west_max-=360;

			if (north_min<min_lat)
				min_lat=north_min;

			if (north_max>max_lat)
				max_lat=north_max;

			if (LonDiff(west_min,min_lon)<0.0)
				min_lon=west_min;

			if (LonDiff(west_max,max_lon)>=0.0)
				max_lon=west_max;
		}

		/* Load any additional SDF files, if required */

		LoadTopoData(max_lon, min_lon, max_lat, min_lat);
	}
 
 // UDT clutter
 LoadUDT   (udt_file);
     


	if (area_mode && topomap==0)
	{

	PlotLRMap(tx_site[0],altitudeLR,ano_filename);

	}

      
	if (map || topomap)
	{


			if (LR.erp==0.0)
				DoPathLoss(mapfile,geo,kml,ngs,tx_site,txsites);
			else
				if (dbm)
					DoRxdPwr(mapfile,geo,kml,ngs,tx_site,txsites);
				else
					DoSigStr(mapfile,geo,kml,ngs,tx_site,txsites);
         }
        
        fprintf(stdout,"|%.5f",north);
        fprintf(stdout,"|%.5f",east);
        fprintf(stdout,"|%.5f",south);
        fprintf(stdout,"|%.5f|",west);
      

        fflush(stdout);

	printf("\n");

	return 0;
}
