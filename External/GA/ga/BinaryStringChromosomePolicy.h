//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef BINARYSTRINGCHROMOSOMEPOLICY_H
#define	BINARYSTRINGCHROMOSOMEPOLICY_H

#include "../util/export.h"
#include "../util/RandomUtil.h"

namespace GapeGa {

class GA_EXPORT BinaryStringChromosomePolicy {
public:
    BinaryStringChromosomePolicy(GarethUtil::RandomUtil & rng_);
    virtual ~BinaryStringChromosomePolicy();
    
    bool mutate(int pos, bool currentValue) const;
    bool initialize(int pos) const;
    bool isAllowSwitch() {return false;}
private:
    GarethUtil::RandomUtil & rng;
    BinaryStringChromosomePolicy(const BinaryStringChromosomePolicy& orig);
    BinaryStringChromosomePolicy & operator=(const BinaryStringChromosomePolicy & other);
};

}


#endif	/* BINARYSTRINGCHROMOSOMEPOLICY_H */

