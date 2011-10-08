
%{
#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include "DiversityPick.h"
%}

%template(EBV_Vect) std::vector< ExplicitBitVect >;

%include "DiversityPick.h";

