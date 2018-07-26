#ifndef __RD_TAUTOMER_CATALOG_UTILS_H__
#define __RD_TAUTOMER_CATALOG_UTILS_H__

#include <GraphMol/RDKitBase.h>
#include "TautomerCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Bond.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TautomerCatalogParams;

// typedef enum {
//	SINGLE,
//	DOUBLE,
//	TRIPLE,
//	AROMATIC,
//} BondType;

// typedef std::vector<ROMol*, std::string, std::string> tautomerTransform;
struct TautomerTransform {
  ROMol* Mol;
  std::vector<Bond::BondType> BondTypes;
  std::vector<int> Charges;
  TautomerTransform(ROMol* mol, std::vector<Bond::BondType> bondtypes,
                    std::vector<int> charges)
      : Mol(mol), BondTypes(bondtypes), Charges(charges) {}
};

std::vector<Bond::BondType> stringToBondType(std::string bond_str);
std::vector<int> stringToCharge(std::string charge_str);

std::vector<TautomerTransform> readTautomers(std::string fileName);
std::vector<TautomerTransform> readTautomers(std::istream& inStream,
                                             int nToRead = -1);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
