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

%{
#include <Numerics/SquareMatrix.h>
#include <Geometry/Transform.h>
#include <Geometry/Transform2D.h>
#include <Geometry/Transform3D.h>
#include <Numerics/SymmMatrix.h>
%}


%include <Numerics/Vector.h>
%include <Numerics/Matrix.h>
%include <Numerics/SquareMatrix.h>
%include <Numerics/SymmMatrix.h>

/* These template definitions need to come after the C++ template definition code, but before the
   instantiation by other modules. */

%template(IntMatrix) RDNumeric::Matrix<int>;
%template(DoubleMatrix) RDNumeric::Matrix<double>;
%template(IntSymmMatrix) RDNumeric::SymmMatrix<int>;
%template(DoubleSymmMatrix) RDNumeric::SymmMatrix<double>;
%template(DoubleSquareMatrix) RDNumeric::SquareMatrix<double>;

%include <Geometry/Transform.h>
%include <Geometry/Transform2D.h>
%include <Geometry/Transform3D.h>

