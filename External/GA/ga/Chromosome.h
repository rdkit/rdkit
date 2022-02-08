//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <string>
#include "../util/export.h"

namespace GapeGa {

/*
 * Base class for all GA chromosomes
 *
 */
class GA_EXPORT Chromosome {
private:
	static int idCounter;
	int chromosomeId;
	double fitness;
	Chromosome (const Chromosome & other) = delete;
	Chromosome & operator = (const Chromosome & other) = delete;

public:
	Chromosome();
	virtual ~Chromosome();

	 virtual Chromosome   & create()=0;

	/*
	 * Rebuilds the chromosome. This ensures that all data structures are
	 * consistent with this chromosome. For example there may be molecular
	 * co-ordinates in a shared molecule object that need to be updated.
	 */
	virtual double rebuild(Chromosome &c)=0;

	/*
	 * Return true if two chromosmes are equal at the Genetic level.
	 */
	virtual bool equals(const Chromosome &c) const=0;

	/*
	 * Return a distance measure for two chromosomes
	 */
	virtual double distance(const Chromosome &c) const=0;

	/*
	 * Return true if the two chromosomes share a niche.
	 */
	virtual bool sameNiche(const Chromosome &c) const =0;

	/*
	 * Return false if the chromosome is invalid. For example the chromosome
	 * cannot be decoded.
	 */
	virtual bool ok() const =0;

	/*
	 * The genetic information in Chromosome c is copied to this chromosome.
	 */
	virtual void copyGene(Chromosome &c) const=0;

	/*
	 * Return a summary of the fitness.
	 */
	virtual std::string fitnessInfo() const = 0;

	/*
	 * Return a summary of the gene.
	 */
	virtual std::string geneInfo() const = 0;

	/*
	 * Determine the fitness of the chromosome
	 */
	virtual void calculateFitness() = 0;

	double getFitness() const {
		return fitness;
	}

	void setFitness(double fitness) {
		this->fitness = fitness;
	}
};

}

#endif /* CHROMOSOME_H_ */
