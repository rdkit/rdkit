//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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

	ROMol *getMetalNof(); // {return metal_nof;}
	ROMol *getMetalNon(); // {return metal_non;}
	void setMetalNof(const ROMol &mol);
	void setMetalNon(const ROMol &mol);

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
