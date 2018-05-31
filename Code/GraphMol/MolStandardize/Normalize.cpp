#include "Normalize.h"
#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

using namespace std;
using namespace RDKit;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

int MAX_RESTARTS = 200;

///////////////////////////////////////////////
const Normalizations_t NORMALIZATIONS[] = {
};
///////////////////////////////////////////////

unsigned int getNumEntries() {
//	static_cast<unsigned int>(sizeof(data) / sizeof(Normalizations_t));
}

Normalization::Normalization( const string &var, const string &transform) {
}; 

ChemicalReaction* Normalization::transform(const string &transform) {
}

RWMol* Normalizer::normalize(const ROMol &mol) {
}

RWMol* Normalizer::normalizeFragment(const ROMol &mol) {
}

RWMol* Normalizer::applyTransform(const ROMol &mol, ChemicalReaction &rule) {
}

} // namespace MolStandardize
} // namespace RDKit

