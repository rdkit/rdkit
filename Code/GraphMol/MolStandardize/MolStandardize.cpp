#include "MolStandardize.h"
//#include <GraphMol/MolStandardize/MolStandardize.h>
////#include <GraphMol/RDKitBase.h>
#include <iostream>
//#include <GraphMol/ROMol.h>
using namespace std;

namespace RDKit{

CleanupParameters::CleanupParameters(){}

namespace MolStandardize{
bool cleanup(RWMol &mol, const CleanupParameters &params){	
	return true;
}
} // end of namespace MolStandardize
} // end of namespace RDKit

