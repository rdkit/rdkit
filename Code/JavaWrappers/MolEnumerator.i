/*
 *
 *  Copyright (c) 2020-2021, Greg Landrum and T5 Informatics GmbH
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


%ignore RDKit::MolEnumerator::detail::idxPropName;
%ignore RDKit::MolEnumerator::detail::preserveOrigIndices;
%ignore RDKit::MolEnumerator::detail::removeOrigIndices;

%ignore RDKit::MolEnumerator::MolEnumeratorOp;
%ignore RDKit::MolEnumerator::PositionVariationOp;
%ignore RDKit::MolEnumerator::LinkNodeOp;
%ignore RDKit::MolEnumerator::RepeatUnitOp;

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
RDKit::MolEnumerator::MolEnumeratorParams getRepeatUnitParams(){
    RDKit::MolEnumerator::MolEnumeratorParams res;
    res.dp_operation.reset(new RDKit::MolEnumerator::RepeatUnitOp());
    return res;
}

%}
%include<GraphMol/MolEnumerator/MolEnumerator.h>
