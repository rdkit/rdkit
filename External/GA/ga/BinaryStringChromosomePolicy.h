/** 
 * File:   BinaryStringChromosomePolicy.h
 * Author: Gareth Jones
 *
 * Created on April 11, 2013, 10:24 PM
 */


#ifndef BINARYSTRINGCHROMOSOMEPOLICY_H
#define	BINARYSTRINGCHROMOSOMEPOLICY_H

#include "StringChromosomeBase.h"
#include "../util/RandomUtil.h"

namespace GapeGa {

class BinaryStringChromosomePolicy {
public:
    BinaryStringChromosomePolicy(GarethUtil::RandomUtil & rng_);
    virtual ~BinaryStringChromosomePolicy();
    
    bool mutate(int pos, bool currentValue) const;
    bool initialize(int pos) const;
    bool isAllowSwitch() {return false;}
private:
    GarethUtil::RandomUtil & rng;
    BinaryStringChromosomePolicy(const BinaryStringChromosomePolicy& orig);
    BinaryStringChromosomePolicy & operator =(const BinaryStringChromosomePolicy & other);
};

}


#endif	/* BINARYSTRINGCHROMOSOMEPOLICY_H */

