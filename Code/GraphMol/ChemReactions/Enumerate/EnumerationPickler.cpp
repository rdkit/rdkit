//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "EnumerationPickler.h"
#include "CartesianProduct.h"
#include "RandomSample.h"
#include "RandomSampleAllBBs.h"

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {

std::string GetClass(const EnumerationStrategyBase *en) {
  if (dynamic_cast<const CartesianProductStrategy *>(en)) {
    return "-->cartesian";
  }
  if (dynamic_cast<const RandomSampleStrategy *>(en)) {
    return "-->random";
  }
  if (dynamic_cast<const RandomSampleAllBBsStrategy *>(en)) {
    return "-->randombbs";
  }
  return "Unknown!";
}

namespace EnumerationStrategyPickler {

void pickle(const boost::shared_ptr<EnumerationStrategyBase> &enumerator,
            std::ostream &ss) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_oarchive ar(ss);
  ar & enumerator;
#else
  RDUNUSED_PARAM(enumerator);
  RDUNUSED_PARAM(ss);
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

void pickle(const boost::shared_ptr<EnumerationStrategyBase> &enumerator,
            std::string &s) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  std::stringstream ss;
  pickle(enumerator, ss);
  s = ss.str();
#else
  RDUNUSED_PARAM(enumerator);
  RDUNUSED_PARAM(s);
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

boost::shared_ptr<EnumerationStrategyBase> fromPickle(std::istream &pickle) {
  boost::shared_ptr<EnumerationStrategyBase> enumerator;
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_iarchive ar(pickle);
  ar & enumerator;
  return enumerator;
#else
  RDUNUSED_PARAM(pickle);
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

boost::shared_ptr<EnumerationStrategyBase> fromPickle(
    const std::string &pickle) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  std::stringstream ss(pickle);
  return fromPickle(ss);
#else
  RDUNUSED_PARAM(pickle);
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
  return boost::shared_ptr<EnumerationStrategyBase>();
#endif
}
}  // namespace EnumerationStrategyPickler
}  // namespace RDKit
