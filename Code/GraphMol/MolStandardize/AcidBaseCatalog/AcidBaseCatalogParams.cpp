#include "AcidBaseCatalogParams.h"
#include "AcidBaseCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>

namespace RDKit {
namespace MolStandardize {

AcidBaseCatalogParams::AcidBaseCatalogParams(const std::string &acidBaseFile) {
	d_pairs.clear();
	d_pairs = readPairs(acidBaseFile);
}

AcidBaseCatalogParams::AcidBaseCatalogParams(const AcidBaseCatalogParams &other) {
	d_typeStr = other.d_typeStr;
	d_pairs.clear();

	const std::vector<std::pair<ROMol*, ROMol*>> &abpairs = other.getPairs();
	for (auto &pairi : abpairs) {
		d_pairs.push_back( std::pair<ROMol*, ROMol*>(pairi.first, pairi.second) );
	}
}

AcidBaseCatalogParams::~AcidBaseCatalogParams() {}

const std::vector<std::pair<ROMol*, ROMol*>> &AcidBaseCatalogParams::getPairs() const {
  return d_pairs;
}

const std::pair<ROMol*, ROMol*> AcidBaseCatalogParams::getPair(unsigned int fid) const {
  URANGE_CHECK(fid, d_pairs.size());
  // return d_pairs[fid];
  return d_pairs[fid];//.get();
}

void AcidBaseCatalogParams::toStream(std::ostream &ss) const {
	ss << d_pairs.size() << "\n";
}

std::string AcidBaseCatalogParams::Serialize() const {
	std::stringstream ss;
	toStream(ss);
	return ss.str();
}

void AcidBaseCatalogParams::initFromStream(std::istream &ss) {}

void AcidBaseCatalogParams::initFromString(const std::string &text) {}

} // namespace MolStandardize
} // namespace RDKit
