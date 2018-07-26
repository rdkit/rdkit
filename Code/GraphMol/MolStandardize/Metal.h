#ifndef __RD_METAL_H__
#define __RD_METAL_H__

#include <GraphMol/ROMol.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {
class MetalDisconnector {
 public:
  MetalDisconnector();
  MetalDisconnector(const MetalDisconnector &other);
  ~MetalDisconnector();

  // TODO inputting ROMol doesn't work
  ROMol *disconnect(const ROMol &mol);
  // overload
  // modifies the molecule in place
  void disconnect(RWMol &mol);  // static?
 private:
  ROMOL_SPTR metal_nof;
  ROMOL_SPTR metal_non;

};  // class Metal
}  // namespace MolStandardize
}  // namespace RDKit
#endif
