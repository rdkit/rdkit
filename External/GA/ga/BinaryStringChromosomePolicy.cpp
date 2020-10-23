/* 
 * File:   BinaryStringChromosomePolicy.cpp
 * Author: gjones
 * 
 * Created on April 11, 2013, 10:24 PM
 */

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
