/*****************************************************************************
*  COST231-HATA MODEL for Signal Server by Alex Farrant                         *  
*  30 December 2013										                     *
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

using namespace std;

double CostHataLinkdB(float f,float h_B, float h_M, float d){
/*
COST HATA URBAN model
Frequency 1500 to 2000MHz
h_B = Base station height 30 to 200m
h_M = Mobile station height 1 to 10m
Distance 1-20km
*/

	int C = 0; // 0dB for suburban
	
	float lh_M = log10(11.75*h_M);
	float C_H = 3.2*lh_M*lh_M-4.97;
	
	float logf = log10(f);
	
    double dbloss = 46.3 + (33.9 * logf) - (13.82 * log10(h_B)) - C_H + (44.9 - 6.55 * log10(h_B)) * log10(d) + C;  


	return dbloss;
}
