//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "BinaryStringChromosomePolicy.h"

namespace GapeGa {

BinaryStringChromosomePolicy::BinaryStringChromosomePolicy(GarethUtil::RandomUtil & rng_) : rng(rng_){
}

BinaryStringChromosomePolicy::~BinaryStringChromosomePolicy() {
}

bool BinaryStringChromosomePolicy::mutate(int pos, bool currentValue) const {
    return ! currentValue;
}
    
bool BinaryStringChromosomePolicy::initialize(int pos) const {
    return rng.randomBoolean();
}

}
