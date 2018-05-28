#ifndef __RD_MOLSTANDARDIZE_H__
#define __RD_MOLSTANDARDIZE_H__

namespace RDKit{
class RWMol;
class ROMol;

struct CleanupParameters{
	// TODO
	// std::vector<std::string> normalisations;
	// std::vector<std::string> acidBasePairs;
	// std::vector<std::string> chargeCorrections;
	// std::vector<std::string> tautomerTransforms;
	// std::vector<std::string> TautomerScores;
	int maxRestarts; // The maximum number of times to attempt to apply the series of normalizations (default 200).
	int maxTautomers; // The maximum number of tautomers to enumerate (default 1000).
	bool preferOrganic; // Whether to prioritize organic fragments when choosing fragment parent (default False).

	CleanupParameters()
		: // TODO 
  		  // normalisations(NORMALISATIONS)
		  // acidBasePairs()
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

void stereoParent(RWMol &mol, const CleanupParameters &params);

void isotopeParent(RWMol &mol, const CleanupParameters &params);

void chargeParent(RWMol &mol, const CleanupParameters &params);

void superParent(RWMol &mol, const CleanupParameters &params);

}; // MolStandardize
}
#endif
