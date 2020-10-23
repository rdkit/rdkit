/*
 * BinaryTestGa.h
 *
 * A test GA to solve the binary f6 function documented in Davis's book.
 *
 *  Created on: Apr 26, 2013
 *      Author: gjones
 */

#ifndef BINARYTESTGA_H_
#define BINARYTESTGA_H_

#include <vector>
#include <memory>
#include "BinaryTestGaChromosome.h"
#include "GaBase.h"
#include "GaOperation.h"
#include "LinkedPopLinearSel.h"
#include "../util/Reporter.h"

namespace GapeGa {

using namespace std;

class BinaryTestGa;
// typdef for the population
typedef LinkedPopLinearSel<BinaryTestGaChromosome, BinaryTestGa> BinaryTestGaPopulation;

class BinaryTestGa: public GaBase {
private:
	BinaryTestGa(const BinaryTestGa & other) = delete;
	BinaryTestGa & operator =(const BinaryTestGa & other) = delete;

	static void mutate(const vector<shared_ptr<BinaryTestGaChromosome> > & parents,
			vector<shared_ptr<BinaryTestGaChromosome> > & children);
	static void crossover(const vector<shared_ptr<BinaryTestGaChromosome> > & parents,
			vector<shared_ptr<BinaryTestGaChromosome> > & children);

	vector<shared_ptr<GaOperation<BinaryTestGaChromosome> > > operations;
	void createOperations();

	BinaryStringChromosomePolicy binaryStringChromosomePolicy;
	// need to have population as a pointer as we want to change popsize etc before creating population
	unique_ptr<BinaryTestGaPopulation> population;
public:
	BinaryTestGa() : binaryStringChromosomePolicy(getRng()) {
		setPopsize(20);
		createOperations();
		population = unique_ptr<BinaryTestGaPopulation>(new BinaryTestGaPopulation(*this));
	}

	virtual ~BinaryTestGa() {
	}

	const vector<shared_ptr<GaOperation<BinaryTestGaChromosome> > > & getOperations() const {
		return operations;
	}

	shared_ptr<BinaryTestGaChromosome> createChromosome();

	void run();

};

}

#endif /* BINARYTESTGA_H_ */
