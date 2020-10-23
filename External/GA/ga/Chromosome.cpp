/*
 * Chromosome.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: Gareth Jones
 */

#include "Chromosome.h"

namespace GapeGa {

int Chromosome::idCounter = 0;

Chromosome::Chromosome() :  chromosomeId(idCounter++), fitness(0){
}

Chromosome::~Chromosome() {
}

}
