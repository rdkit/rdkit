/*
 *
 *  Copyright (c) 2020, Greg Landrum and T5 Informatics GmbH
 *  All rights reserved.
 *
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
 *
 */
%{
#include <GraphMol/MolEnumerator/MolEnumerator.h>
%}


#if 0
%ignore RDKit::MolEnumerator::MolEnumeratorOp::operator();
%ignore RDKit::MolEnumerator::MolEnumeratorOp::copy;
%ignore RDKit::MolEnumerator::PositionVariationOp::operator();
%ignore RDKit::MolEnumerator::PositionVariationOp::copy;
%ignore RDKit::MolEnumerator::LinkNodeOp::operator();
%ignore RDKit::MolEnumerator::LinkNodeOp::copy;
#endif

%ignore RDKit::MolEnumerator::MolEnumeratorOp;
%ignore RDKit::MolEnumerator::PositionVariationOp;
%ignore RDKit::MolEnumerator::LinkNodeOp;

%inline %{

RDKit::MolEnumerator::MolEnumeratorParams getLinkNodeParams(){
    RDKit::MolEnumerator::MolEnumeratorParams res;
    res.dp_operation.reset(new RDKit::MolEnumerator::LinkNodeOp());
    return res;
}
RDKit::MolEnumerator::MolEnumeratorParams getPositionVariationParams(){
    RDKit::MolEnumerator::MolEnumeratorParams res;
    res.dp_operation.reset(new RDKit::MolEnumerator::PositionVariationOp());
    return res;
}

%}
%include<GraphMol/MolEnumerator/MolEnumerator.h>
