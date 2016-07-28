//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include "../RDKitBase.h"

namespace RDKit {
    namespace StructureCheck {
//TODO: ring_list type
        void createRingList(const ROMol &mol, std::vector<unsigned> &ring_list);
        void CombineRings(const ROMol &mol, std::vector<unsigned> &ring_list);
    }
}

