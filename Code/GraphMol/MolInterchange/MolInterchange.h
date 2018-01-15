//
//  Copyright (C) 2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_MOLINTERCHANGE_H_JAN2018
#define RD_MOLINTERCHANGE_H_JAN2018

#include <string>
#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace RDKit {

class RWMol;

namespace MolInterchange {
// \brief construct molecules from MolJSON data in a stream
/*!
 *   \param inStream - stream containing the data *
 */
std::vector<boost::shared_ptr<RWMol>> JSONDataStreamToMols(
    std::istream *inStream);

// \brief construct a molecule from an MDL mol block
/*!
 *   \param jsonBlock - string containing the mol block
 */
std::vector<boost::shared_ptr<RWMol>> JSONDataToMols(
    const std::string &jsonBlock);

std::string MolsToJSONData(const std::vector<const ROMol *> &mols);
std::string MolToJSONData(const ROMol &mol) {
  std::vector<const ROMol *> ms{&mol};
  return MolsToJSONData(ms);
};

}  // end of namespace MolInterchange
}  // end of namespace RDKit

#endif
