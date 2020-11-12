//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef BINARYSTRINGCHROMOSOME_H_
#define BINARYSTRINGCHROMOSOME_H_

#include "StringChromosomeBase.h"
#include "BinaryStringChromosomePolicy.h"
#include "IntegerStringChromosomePolicy.h"

namespace GapeGa {

template<typename T, typename ChromosomePolicy>
class StringChromosome: public StringChromosomeBase<T, ChromosomePolicy> {
private:
	// we should never use the default implementation
	StringChromosome() = delete;
};

template<> class StringChromosome<bool, BinaryStringChromosomePolicy> : public StringChromosomeBase<
		bool, BinaryStringChromosomePolicy> {
private:
	StringChromosome(const StringChromosome & other) = delete;
	StringChromosome & operator =(const StringChromosome & other) = delete;
public:
	StringChromosome(int length_, GarethUtil::RandomUtil & rng_,
			BinaryStringChromosomePolicy & chromosomePolicy_) :
			StringChromosomeBase(length_, rng_, chromosomePolicy_) {
		;
	}

	int decodeToInt(int start, int nBits) const;
	int decodeByte(int byteNo) const {
		return decodeToInt(byteNo*8, 8);
	}

	void test();
};



using BinaryStringChromosome = StringChromosome<bool, BinaryStringChromosomePolicy> ;
using IntegerStringChromosome = StringChromosomeBase<int, IntegerStringChromosomePolicy>;

}

#endif /* BINARYSTRINGCHROMOSOME_H_ */
