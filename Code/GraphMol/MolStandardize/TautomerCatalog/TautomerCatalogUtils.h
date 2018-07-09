#ifndef __RD_TAUTOMER_CATALOG_UTILS_H__
#define __RD_TAUTOMER_CATALOG_UTILS_H__

#include <GraphMol/RDKitBase.h>
#include "TautomerCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TautomerCatalogParams;

//typedef std::vector<ROMol*, std::string, std::string> tautomerTransform;
struct TautomerTransform {
				ROMol* Mol;
				std::string Bonds;
				std::string Charges;
				TautomerTransform(ROMol* mol, std::string bonds,
												std::string charges)
								: Mol(mol), Bonds(bonds), Charges(charges) {}
};

std::vector<TautomerTransform> readTautomers(std::string fileName);
std::vector<TautomerTransform> readTautomers(std::istream &inStream, int nToRead = -1);

} // namespace MolStandardize
} // namespace RDKit

#endif
