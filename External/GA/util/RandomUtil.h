/*
 * Utils.h
 *
 *  Created on: Mar 28, 2013
 *      Author: Gareth Jones
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <cstdint>
#include <random>

namespace GarethUtil {

/*
 * Singleton class to hold a random number generator
 *
 */
class RandomUtil {
public:

	RandomUtil(const RandomUtil & rhs) = delete;
	RandomUtil & operator =(const RandomUtil & rhs) = delete;
	RandomUtil(RandomUtil && rhs) = delete;
	RandomUtil & operator =(RandomUtil && rhs) = delete;

	/*
	 * return a random number between 0 and 1
	 */
	double normalRand();

	/*
	 * return a random integer between top and bottom
	 */
	int randomInt(int bottom, int top);

	bool randomBoolean();

	void seed(uint32_t seed);

	/**
	 * Get singleton
	 */
	static RandomUtil &getInstance();
private:
	RandomUtil();
	virtual ~RandomUtil();

	std::mt19937 rng;
	std::uniform_real_distribution<double> realDistribution;

};

}
#endif /* UTILS_H_ */
