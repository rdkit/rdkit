/*
 * StringChromosome.cpp
 *
 *  Created on: May 4, 2013
 *      Author: gjones
 */

#include "StringChromosome.h"

namespace GapeGa {

/**
 * A decode to integer function for binary string chromosomes
 *
 * @param start
 * @param nBits
 * @return
 */
int StringChromosome<bool, BinaryStringChromosomePolicy>::decodeToInt(
		int start, int nBits) const {

	assert(length >= (start+nBits));

	int mask = 1, result = 0;
	bool *ptr = string.get()+start;
	for (int i =0; i < nBits; i++, mask <<= 1) {
		if (*ptr++ == true)
			result |= mask;
	}
	return result;
}

}
