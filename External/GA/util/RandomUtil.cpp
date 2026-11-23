//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RandomUtil.h"
#include <mutex>

using std::mt19937;
using std::uniform_real_distribution;

namespace GarethUtil {

namespace {
// Parallel GA runs share the singleton engine and distribution. Use the same
// mutex for seeding and draws.
std::mutex randomMutex;
}

// RandomUtil::RandomUtil(uint32_t seed_) :
//		realDistribution(0, 1){
//	seed(seed_);
//}

RandomUtil::RandomUtil() : realDistribution(0, 1) {}

RandomUtil::~RandomUtil() {}

void RandomUtil::seed(uint32_t seed_) {
  const std::lock_guard<std::mutex> lock(randomMutex);
  rng.seed(seed_);
}

double RandomUtil::normalRand() {
  const std::lock_guard<std::mutex> lock(randomMutex);
  return realDistribution(rng);
}

int RandomUtil::randomInt(int bottom, int top) {
  const std::lock_guard<std::mutex> lock(randomMutex);
  int step = rng() % (top - bottom);
  return bottom + step;
}

bool RandomUtil::randomBoolean() { return normalRand() < 0.5; }

RandomUtil &RandomUtil::getInstance() {
  static RandomUtil rng;
  return rng;
}

}  // namespace GarethUtil
