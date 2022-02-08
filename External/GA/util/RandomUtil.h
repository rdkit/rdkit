//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef UTILS_H_
#define UTILS_H_

#include <cstdint>
#include <random>
#include "export.h"

namespace GarethUtil {

/*
 * Singleton class to hold a random number generator
 *
 */
class GA_EXPORT RandomUtil {
 public:
  RandomUtil(const RandomUtil &rhs) = delete;
  RandomUtil &operator=(const RandomUtil &rhs) = delete;
  RandomUtil(RandomUtil &&rhs) = delete;
  RandomUtil &operator=(RandomUtil &&rhs) = delete;

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

}  // namespace GarethUtil
#endif /* UTILS_H_ */
