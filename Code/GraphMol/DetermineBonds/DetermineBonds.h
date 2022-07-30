//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#pragma once
#include <GraphMol/RDKitBase.h>

namespace RDKit {

void determineConnectivity(RWMol &mol, bool useHuckel=false, int charge=0, double covFactor=1.3);

void connectivityHuckel(RWMol &mol, int charge=0);

void connectivityVdW(RWMol &mol, double covFactor=1.3);

}
