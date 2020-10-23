/*
 * IntegerStringChromosomePolicy.h
 *
 *  Created on: Apr 12, 2013
 *      Author: gjones
 */

#ifndef INTEGERSTRINGCHROMOSOMEPOLICY_H_
#define INTEGERSTRINGCHROMOSOMEPOLICY_H_

#include "../util/RandomUtil.h"
#include "StringChromosomeBase.h"

namespace GapeGa {

class IntegerStringChromosomePolicy {
public:
	IntegerStringChromosomePolicy(GarethUtil::RandomUtil & rng_, int s);
	virtual ~IntegerStringChromosomePolicy();

	int mutate(int pos, int currentValue) const;
	int initialize(int pos) const;
	bool isAllowSwitch() {
		return false;
	}

	void setMax(int max);
	void setMax(int pos, int max);
	void setAllowNulls(bool allow);
	void setAllowNulls(int pos, bool allow);

	int getSize() const {
		return size;
	}

private:
	GarethUtil::RandomUtil & rng;
	const int size;
	IntegerStringChromosomePolicy(const IntegerStringChromosomePolicy& orig) = delete;
	IntegerStringChromosomePolicy & operator =(
			const IntegerStringChromosomePolicy & other) = delete;

	int * const maxs;
	bool * const allowNulls;
};

}

#endif /* INTEGERSTRINGCHROMOSOMEPOLICY_H_ */
