/*
 * BinaryTestGaChromosome.cpp
 *
 *  Created on: May 4, 2013
 *      Author: gjones
 */

#include "BinaryTestGaChromosome.h"

#include <boost/format.hpp>
#include <cassert>
#include <cmath>
#include <string>

#include "StringChromosomeBase.h"

namespace GapeGa {

using namespace std;

/**
 * binary f6 function (as decscribed in Davis, Handbook of Genetic Algorithms, pg 8).
 * @return
 */
double BinaryTestGaChromosome::score() {
	static double scale = 200.0 / (pow(2.0, 22) - 1);

	int xint = decodeToInt(0, 22);
	int yint = decodeToInt(22, 22);
	xVal = static_cast<double>(xint);
	yVal = static_cast<double>(yint);
	xVal *= scale;
	yVal *= scale;
	assert(xVal <= 200 && xVal >= 0);
	assert(yVal <= 200 && yVal >= 0);
	xVal -= 100;
	yVal -= 100;
	double sqr = xVal * xVal + yVal * yVal;
	double s = sin(sqrt(sqr));
	double num = s * s - 0.5;
	double denom = 1.0 + 0.0001 * sqr;
	denom = denom * denom;
	fitness = 0.5 - num / denom;
	return fitness;
}

string BinaryTestGaChromosome::info() const {
	boost::format format = boost::format(
				"Fit %6.3f x %6.3f y %6.3f : ") % fitness % xVal % yVal;
	return format.str()+ geneInfo();
}

}
