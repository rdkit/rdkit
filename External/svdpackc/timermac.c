/*************************************************************************
                           (c) Copyright 1993
                        University of Tennessee
                          All Rights Reserved                          
 *************************************************************************/
/*--------------------------------------------------------------------------
 *  Timing Sample Program
 * 
 *	This programs uses the "TickCount" Macintosh Tool which gives the 
 *  system time in 1/60th of a second. It computes the elapsed time 
 *  and converts it into min/sec units.
 *
 *------------------------------------------------------------------------*/
 
#include	<stdio.h> 
#include	<Packages.h>

float timer()
{
   return( (float) TickCount() / 60.0);
 } 
