#ifndef __RD_MOLSTANDARDIZE_H__
#define __RD_MOLSTANDARDIZE_H__

namespace RDKit{
class RWMol;
class ROMol;

struct CleanupParameters{
	CleanupParameters();
};

namespace MolStandardize{
	
bool cleanup(RWMol &mol, const CleanupParameters &params);

}; // MolStandardize
}
#endif
