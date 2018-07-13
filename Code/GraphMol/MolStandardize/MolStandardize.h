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
	std::string tautomerTransforms;
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
		  tautomerTransforms(rdbase + "/Code/GraphMol/MolStandardize/TautomerCatalog/test_data/tautomerTransforms.in"),
		  // TautomerScores()
		  maxRestarts(200),
		  maxTautomers(1000),
		  preferOrganic(false)
	{}

};

namespace MolStandardize{
	
RWMol* cleanup(const RWMol &mol, const CleanupParameters &params);

void tautomerParent(RWMol &mol, const CleanupParameters &params);

RWMol* fragmentParent(const RWMol &mol, const CleanupParameters &params, 
		bool skip_standardize = false);

void stereoParent(RWMol &mol, const CleanupParameters &params);

void isotopeParent(RWMol &mol, const CleanupParameters &params);

RWMol* chargeParent(const RWMol &mol, const CleanupParameters &params, 
								bool skip_standardize=false);

void superParent(RWMol &mol, const CleanupParameters &params);

RWMol* normalize(const RWMol *mol, const CleanupParameters &params);

RWMol* reionize(const RWMol *mol, const CleanupParameters &params);

std::string standardizeSmiles(const std::string &smiles);

std::vector<std::string> enumerateTautomerSmiles(const std::string &smiles, 
								const CleanupParameters &params);
}; // MolStandardize
}
#endif
