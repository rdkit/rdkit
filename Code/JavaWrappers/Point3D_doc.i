/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%typemap(javaimports) RDGeom::Point3D "
/** Class representing 3D point with x, y and z coordinates. */"

%javamethodmodifiers RDGeom::Point3D::angleTo 	( 	const Point3D &  	other 	 )  	const "
/**
<p>
determines the angle between a vector to this point from the origin and a vector to the other point.
<p>
The angle is unsigned: the results of this call will always be between 0 and M_PI
*/
public";

%javamethodmodifiers RDGeom::Point3D::crossProduct 	( 	const Point3D &  	other 	 )  	const "
/**
<p>
Cross product of this point with the another point.
<p>
The order is important here The result is 'this' cross with 'other' not (other x this)
*/
public";

%javamethodmodifiers RDGeom::Point3D::directionVector 	( 	const Point3D &  	other 	 )  	const "
/**
<p>
<p>
@return
a normalized direction vector from this point to another.
.
*/
public";

%javamethodmodifiers RDGeom::Point3D::getPerpendicular 	( 		 )  	const "
/**
<p>
Get a unit perpendicular from this point (treating it as a vector):
*/
public";

%javamethodmodifiers RDGeom::Point3D::signedAngleTo 	( 	const Point3D &  	other 	 )  	const "
/**
<p>
determines the signed angle between a vector to this point from the origin and a vector to the other point.
<p>
The results of this call will be between 0 and M_2_PI
<p>

*/
public";

