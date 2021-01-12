//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "IntegerStringChromosomePolicy.h"

namespace GapeGa {

using namespace GarethUtil;

IntegerStringChromosomePolicy::IntegerStringChromosomePolicy(RandomUtil & rng_,
		int size_) :
		rng(rng_), size { size_ }, maxs { new int[size] }, allowNulls {
				new bool[size] } {
					setMax(10);
					setAllowNulls(false);
}

IntegerStringChromosomePolicy::~IntegerStringChromosomePolicy() {
	delete[] maxs;
	delete[] allowNulls;
}

int IntegerStringChromosomePolicy::mutate(int pos, int currentValue) const {
	assert(pos >= 0 && pos < size);
	int max = maxs[pos];
	int min = allowNulls[pos] ? -1 : 0;
	int val = rng.randomInt(min, max - 1);
	if (val >= currentValue) {
		val++;
	}
	return val;
}

int IntegerStringChromosomePolicy::initialize(int pos) const {
	assert(pos >= 0 && pos < size);
	int max = maxs[pos];
	int min = allowNulls[pos] ? -1 : 0;
	return rng.randomInt(min, max);
}

void IntegerStringChromosomePolicy::setMax(int max) {
	for (int i = 0; i < size; i++) {
		maxs[i] = max;
	}
}

void IntegerStringChromosomePolicy::setMax(int pos, int max) {
	assert(pos >= 0 && pos < size);
	maxs[pos] = max;
}

void IntegerStringChromosomePolicy::setAllowNulls(bool allow) {
	for (int i = 0; i < size; i++) {
		allowNulls[i] = allow;
	}

}

void IntegerStringChromosomePolicy::setAllowNulls(int pos, bool allow) {
	assert(pos >= 0 && pos < size);
	allowNulls[pos] = allow;
}

}
