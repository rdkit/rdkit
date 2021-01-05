//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef GABASE_H_
#define GABASE_H_

#include <memory>
#include "../util/RandomUtil.h"
#include "../util/export.h"

namespace GapeGa {

using namespace GarethUtil;
/*
 * Interface to support different GAs. Allows other objects to be able to access
 * the main GA class. Some of these objects may not be defined by GAs- in which
 * case an exception should be thrown.
 *
 */
class GA_EXPORT GaBase {
 private:
  std::string fileName;
  GarethUtil::RandomUtil& rng = RandomUtil::getInstance();
  size_t popsize = 100;
  double selectionPressure = 1.1;
  GaBase(const GaBase& other) = delete;
  GaBase& operator=(const GaBase& other) = delete;

 public:
  GaBase(){};
  virtual ~GaBase(){};

  double getSelectionPressure() const { return selectionPressure; }

  size_t getPopsize() const { return popsize; }

  GarethUtil::RandomUtil& getRng() { return rng; }

 protected:
  void setSelectionPressure(double selectionPressure) {
    this->selectionPressure = selectionPressure;
  }

  void setPopsize(size_t popsize) { this->popsize = popsize; }
};

}  // namespace GapeGa

#endif /* GABASE_H_ */
