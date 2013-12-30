/*****************************************************************************
*  HATA MODEL for Signal Server by Alex Farrant                              *
*  30 December 2013								                             *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.								   								     *
* 									    									 *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.													     *
*									     									 */

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

#define PI 3.14159265

/* Acute Angle from Rx point to an obstacle of height (opp) and distance (adj) */ 
double incidenceAngle(double opp, double adj){
return atan2(opp,adj) * 180 / PI;
}

/* 
Knife edge diffraction:
This custom function adds to the overall path loss based upon the obstacle 
angle (from receive point) and the frequency. Loss increases with angle and frequency.
This is not a recognised formula like Huygens, rather it is a 
compromise for increased speed which adds a realistic diffraction effect.
*/
double ked(double freq, double elev[], double rxh, double dkm){
double obh,obd,rxobaoi=0,d;

obh=0;
obd=0;

dkm=dkm*1000; // KM to metres

	for(int n=2;n<(dkm/elev[1]);n++){
	
	d = (n-2)*elev[1]; // no of points * delta = km
	
	    //Find dip(s)
		if(elev[n] > 0 && elev[n]<(obh+10)){
		
		// Angle from Rx point to obstacle
		rxobaoi = incidenceAngle((obh-(elev[n]+rxh)),d-obd);
		} else{
		// Line of sight or higher
		rxobaoi=0;
		}
	
		//note the highest point
		if(elev[n]>obh){
		obh=elev[n];
		obd=d; 
		}
		
	}
	
if(rxobaoi >= 0){
return (rxobaoi * 3) / (300/freq); // Exaggerate diffraction angle and divide by wavelength (m)
}else{
return 0;
}

}

double HataLinkdB(float f,float h_B, float h_M, float d, int mode){
/*
HATA URBAN model for cellular planning
Frequency (MHz) 150 to 1500MHz
Base station height 30-200m
Mobile station height 1-10m
Distance 1-20km

mode 1 = URBAN
mode 2 = SUBURBAN
mode 3 = OPEN
*/

	float lh_M = log10(11.75*h_M);
	float C_H = 3.2*lh_M*lh_M-4.97;
	
	float logf = log10(f);
	
	float L_u = 69.55 + 26.16*logf - 13.82*log10(h_B) - C_H + (44.9 - 6.55*log10(h_B))*log10(d);

	if(!mode || mode==1){
	return L_u; //URBAN
	}
	
	if(mode==2){ //SUBURBAN
	float logf_28 = log10(f/28);
	return L_u - 2*logf_28*logf_28 - 5.4;
	}

	if(mode==3){ //OPEN
	return L_u - 4.78*logf*logf + 18.33*logf - 40.94;
	}

	return 0;
}