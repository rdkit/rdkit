#ifndef __RD_MOLSTANDARDIZE_H__
#define __RD_MOLSTANDARDIZE_H__

#include <string>
#include <GraphMol/RDKitBase.h>

namespace RDKit{
class RWMol;
class ROMol;

struct CleanupParameters{
	// TODO
	std::string rdbase = std::getenv("RDBASE");
	std::string normalizations;
	std::string acidbaseFile;
	// std::vector<std::string> chargeCorrections;
	// std::vector<std::string> tautomerTransforms;
	// std::vector<std::string> TautomerScores;
	int maxRestarts; // The maximum number of times to attempt to apply the series of normalizations (default 200).
	int maxTautomers; // The maximum number of tautomers to enumerate (default 1000).
	bool preferOrganic; // Whether to prioritize organic fragments when choosing fragment parent (default False).

	CleanupParameters()
		: // TODO 
//			normalizations(""),//this->DEFAULT_TRANSFORMS),
			normalizations(rdbase + "/Code/GraphMol/MolStandardize/TransformCatalog/test_data/normalizations.txt"),
		  acidbaseFile(rdbase + "/Code/GraphMol/MolStandardize/AcidBaseCatalog/test_data/acid_base_pairs.txt"),
		  // chargeCorrections()
		  // tautomerTransforms()
		  // TautomerScores()
		  maxRestarts(200),
		  maxTautomers(1000),
		  preferOrganic(false)
	{}

};

namespace MolStandardize{
	
bool cleanup(RWMol &mol, const CleanupParameters &params);

void tautomerParent(RWMol &mol, const CleanupParameters &params);

void fragmentParent(RWMol &mol, const CleanupParameters &params, 
		bool skip_standardize = false);

void stereoParent(RWMol &mol, const CleanupParameters &params);

void isotopeParent(RWMol &mol, const CleanupParameters &params);

void chargeParent(RWMol &mol, const CleanupParameters &params, 
								bool skip_standardize=false);

void superParent(RWMol &mol, const CleanupParameters &params);

void normalize(RWMol &mol, const CleanupParameters &params);

void reionize(RWMol &mol, const CleanupParameters &params);

std::string standardizeSmiles(const std::string &smiles);

}; // MolStandardize
}
#endif
