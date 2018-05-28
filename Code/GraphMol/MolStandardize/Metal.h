#ifndef __RD_METAL_H__
#define __RD_METAL_H__

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{
class MetalDisconnector{
	public:
		MetalDisconnector();
		ROMol* disconnect(const ROMol &mol);
		// overload
		// modifies the molecule in place
		void disconnect(RWMol &mol);
	private:
		ROMol* metal_nof;
		ROMol* metal_non;

}; // class Metal
} // namespace MolStandardize
} // namespace RDKit
#endif
