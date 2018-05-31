#ifndef __RD_NORMALIZE_H__
#define __RD_NORMALIZE_H__

#include <string>
using namespace std;

namespace RDKit{
class RWMol;
class ROMol;
class ChemicalReaction;

namespace MolStandardize{

struct Normalizations_t {
	const char* description;
	const char* smarts;
};

unsigned int getNumEntries();

class Normalization {
	public:
		Normalization( const string &var, const string &transform );
		ChemicalReaction* transform(const string &transform);

}; // Normalization class

class Normalizer {
	public:
		Normalizer();
		RWMol* normalize(const ROMol &mol);
	private:
		RWMol* normalizeFragment(const ROMol &mol);
		RWMol* applyTransform(const ROMol &mol, ChemicalReaction &rule);

}; // Normalizer class
}
}

#endif
