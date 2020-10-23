/*
 * Utils.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: Gareth Jones
 */

#include "RandomUtil.h"

using std::mt19937;
using std::uniform_real_distribution;

namespace GarethUtil {

//RandomUtil::RandomUtil(uint32_t seed_) :
//		realDistribution(0, 1){
//	seed(seed_);
//}

RandomUtil::RandomUtil() :
		realDistribution(0, 1){
}

RandomUtil::~RandomUtil() {
}

void RandomUtil::seed(uint32_t seed_) {
	rng.seed(seed_);
}

double RandomUtil::normalRand() {
	return realDistribution(rng);
}


int RandomUtil::randomInt(int bottom, int top) {
	int step = rng() % (top-bottom);
	return bottom+step;
}

bool RandomUtil::randomBoolean() {
	return normalRand()<0.5;
}

RandomUtil &RandomUtil::getInstance() {
	static RandomUtil rng;
	return rng;
}

}
