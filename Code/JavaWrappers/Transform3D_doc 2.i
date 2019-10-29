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

%typemap(javaimports) RDGeom::Transform3D "
/**  */"

%javamethodmodifiers RDGeom::Transform3D::Transform3D 	( 		 )  	"
/**
<p>
Constructor.
<p>
Initialize to an identity matrix transformation. This is a 4x4 matrix that includes the rotation and translation parts see Foley's 'Introduction to Computer Graphics' for the representation
<p>
Operator *= and = are provided by the parent class square matrix. Operator *= needs some explanation, since the order matters. This transform gets set to the combination other and the current state of this transform If this_old and this_new are the states of this object before and after this function we have this_new(point) = this_old(other(point))

*/
public";

%javamethodmodifiers RDGeom::Transform3D::SetRotation 	( 	double  	angle, 		const Point3D &  	axis	  	) 			"
/**
<p>
set the rotation matrix
<p>
The rotation matrix is set to rotation by th specified angle about an arbitrary axis
*/
public";

%javamethodmodifiers RDGeom::Transform3D::SetRotation 	( 	double  	angle, 		AxisType  	axis	  	) 			"
/**
<p>
set the rotation matrix
<p>
The rotation matrix is set to rotation by th specified angle about the specified axis
<p>

*/
public";

