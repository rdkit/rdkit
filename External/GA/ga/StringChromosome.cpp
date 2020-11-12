//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
