#include "Metal.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>

using namespace std;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

MetalDisconnector::MetalDisconnector()
	: metal_nof(SmartsToMol("[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]")),
	  metal_non(SmartsToMol("[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]")){
	  };

ROMol* MetalDisconnector::disconnect(const ROMol &mol){
	auto *res = new RWMol(mol);
	MetalDisconnector::disconnect(*res);
	return static_cast<ROMol *>(res);
}

void MetalDisconnector::disconnect(RWMol &mol){

}

} // namespace MolStandardize
} // namespace RDKit

