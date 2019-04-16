//
//  Copyright (C) 2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLINTERCHANGE_H_JAN2018
#define RD_MOLINTERCHANGE_H_JAN2018

/*! \file MolInterchange.h

\brief Contains the public API for the convertors to/from the commonchem
interchange format

\b Note that this should be considered beta and that the format and API
will very likely change in future releases.

More information about CommonChem is available here:
https://github.com/mcs07/CommonChem
*/

#include <string>
#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace RDKit {

class RWMol;

namespace MolInterchange {

// \brief parameters controlling parsing of MolJSON
struct RDKIT_MOLINTERCHANGE_EXPORT JSONParseParameters {
  bool setAromaticBonds =
      true; /*! toggles setting the BondType of aromatic bonds to Aromatic */
  bool strictValenceCheck =
      false; /*! toggles doing reasonable valence checks */
  bool parseProperties =
      true; /*! toggles extracting molecular properties from the JSON block */
  bool parseConformers =
      true; /*! toggles extracting conformers from the JSON block */
};
static JSONParseParameters defaultJSONParseParameters;

// \brief construct molecules from MolJSON data in a stream
/*!
 *   \param inStream - stream containing the data
 *   \param params   - parsing options
 */
RDKIT_MOLINTERCHANGE_EXPORT std::vector<boost::shared_ptr<ROMol>>
JSONDataStreamToMols(std::istream *inStream, const JSONParseParameters &params =
                                                 defaultJSONParseParameters);

// \brief construct molecules from MolJSON data
/*!
 *   \param jsonBlock - string containing the mol block
 *   \param params   - parsing options
 */
RDKIT_MOLINTERCHANGE_EXPORT std::vector<boost::shared_ptr<ROMol>>
JSONDataToMols(const std::string &jsonBlock,
               const JSONParseParameters &params = defaultJSONParseParameters);

// \brief returns MolJSON for a set of molecules
/*!
 *   \param mols  - the molecules to work with
 */
template <typename T>
RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData(
    const std::vector<T> &mols);

// \brief returns MolJSON for a molecule
/*!
 *   \param mol   - the molecule to work with
 */
template <typename T>
std::string MolToJSONData(const T &mol) {
  std::vector<const T *> ms{&mol};
  return MolsToJSONData(ms);
};

}  // end of namespace MolInterchange
}  // end of namespace RDKit

#endif
