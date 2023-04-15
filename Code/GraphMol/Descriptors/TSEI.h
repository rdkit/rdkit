//
// Created by Jason Biggs on 4/15/23.
//

#include <RDGeneral/export.h>
#ifndef RDKIT_TSEI_H
#define RDKIT_TSEI_H

#include <vector>
#include <string>

namespace RDKit {
class ROMol;
namespace Descriptors {
const std::string TSEIVersion = "1.0.0";

//! calculates MQN descriptors
/*!
  Chenzhong Cao and Li Liu
  "Topological steric effect index and its application"
  CJ Chem Inf Comput Sci 44(2), 678-687, 2004.


  \param mol        the molecule of interest
  \param force      (optional) calculate the values even if they are cached.

  \return a vector with the TSEIs

*/
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTSEI(const ROMol &mol,
                                                      bool force = false);

}  // end of namespace Descriptors
}  // end of namespace RDKit

#endif  // RDKIT_TSEI_H
