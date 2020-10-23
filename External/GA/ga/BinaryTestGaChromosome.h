/**
 * BinaryTestGaChromosome.h
 *
 * A class to represent a Chromosome used in the binary f6 problem
 *
 *  Created on: May 4, 2013
 *      Author: gjones
 */

#ifndef BINARYTESTGACHROMOSOME_H_
#define BINARYTESTGACHROMOSOME_H_

#include "StringChromosome.h"

namespace GapeGa {

class BinaryTestGaChromosome: public BinaryStringChromosome {
private:
	BinaryTestGaChromosome(const BinaryTestGaChromosome & other) = delete;
	BinaryTestGaChromosome & operator =(const BinaryTestGaChromosome & other) = delete;
	double fitness = .0, xVal = .0, yVal = .0;
public:
	BinaryTestGaChromosome(RandomUtil & rng_,
			BinaryStringChromosomePolicy & chromosomePolicy_) :
			BinaryStringChromosome(44, rng_, chromosomePolicy_) {
	}
	virtual ~BinaryTestGaChromosome() {
	}

	double getFitness() const {
		return fitness;
	}

	bool isOk() {
		return true;
	}

	double score();

	double getXVal() const {
		return xVal;
	}

	double getYVal() const {
		return yVal;
	}

	std::string info() const;
};

}

#endif /* BINARYTESTGACHROMOSOME_H_ */
