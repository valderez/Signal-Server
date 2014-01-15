/*****************************************************************************
*  ITU-R P.525 Free Space Path Loss model for Signal Server by Alex Farrant  *  
*  15 January 2014										                     *
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

double FsplLinkdB(float f, float d){
/*
Free Space Path Loss model
Frequency: Any
Distance: Any 
*/
	//MHz to GHz
	f = f / 1000;

	double dbloss = (20 * log10(d)) + (20 * log10(f)) + 92.45;
   
	return dbloss;
}
