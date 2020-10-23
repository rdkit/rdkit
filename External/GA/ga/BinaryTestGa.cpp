/*
 * BinaryTestGa.cpp
 *
 *  Created on: Apr 26, 2013
 *      Author: gjones
 */

#include "BinaryTestGa.h"

namespace GapeGa {

using namespace std;
using namespace GarethUtil;

/**
 * Create the operations for the binary test GA.  We have a equal probability
 * of mutation and one-point crossover.
 */
void BinaryTestGa::createOperations() {
	shared_ptr<GaOperation<BinaryTestGaChromosome> > mutationOperation(
			new GaOperation<BinaryTestGaChromosome>(1, 1, 50.0,
					&BinaryTestGa::mutate));
	shared_ptr<GaOperation<BinaryTestGaChromosome> > crossoverOperation(
			new GaOperation<BinaryTestGaChromosome>(2, 2, 50.0,
					&BinaryTestGa::crossover));
	operations.reserve(2);
	operations.push_back(mutationOperation);
	operations.push_back(crossoverOperation);
}

/**
 * Applies a single mutation operator
 *
 * @param parents
 * @param children
 */
void BinaryTestGa::mutate(
		const std::vector<shared_ptr<BinaryTestGaChromosome> > & parents,
		std::vector<shared_ptr<BinaryTestGaChromosome> > & children) {

	shared_ptr<BinaryTestGaChromosome> parent = parents[0];
	shared_ptr<BinaryTestGaChromosome> child = children[0];

	child->copyGene(*parent);
	child->mutate();
}

/**
 * Applies one-point crossover
 *
 * @param parents
 * @param children
 */
void BinaryTestGa::crossover(
		const vector<shared_ptr<BinaryTestGaChromosome> > & parents,
		vector<shared_ptr<BinaryTestGaChromosome> > & children) {

	shared_ptr<BinaryTestGaChromosome> parent1 = parents[0];
	shared_ptr<BinaryTestGaChromosome> child1 = children[0];
	shared_ptr<BinaryTestGaChromosome> parent2 = parents[1];
	shared_ptr<BinaryTestGaChromosome> child2 = children[1];

	parent1->onePointCrossover(*parent2, *child1, *child2);

}

/**
 * Create a pointer to a new chromosome.
 *
 * @return
 */
shared_ptr<BinaryTestGaChromosome> BinaryTestGa::createChromosome() {
	shared_ptr<BinaryTestGaChromosome> chromosome(
			new BinaryTestGaChromosome(getRng(), binaryStringChromosomePolicy));
	REPORT(Reporter::DEBUG) << "Created chromosome at address "
			<< chromosome.get();
	return chromosome;
}

/**
 * Run the GA
 */
void BinaryTestGa::run() {
	population->create();
	REPORT(Reporter::INFO) << population->info();
	REPORT(Reporter::DETAIL) << population->populationInfo();

	int nOps = 0;
	while (nOps < 10000) {
		population->iterate();
		nOps++;
		if (nOps % 1000 == 0) {
			REPORT(Reporter::INFO) << population->info();
		}
	}
	const shared_ptr<BinaryTestGaChromosome> best = population->getBest();
	REPORT(Reporter::INFO) << "Best solution " << best->info();
	REPORT(Reporter::DETAIL) << population->populationInfo();

}

}
