//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#ifndef _RD_CHEMTRANSFORMS_RGROUPDECOMP_H__
#define _RD_CHEMTRANSFORMS__RGROUPDECOMPH__

#include <vector>
#include <boost/smart_ptr.hpp>

#include "../RDKitBase.h"

namespace RDKit{ 
    enum RGoupDecompositionSubstituentHandling {
        PASS,   // the substituents through anyway
        IGNORE, // ignore the substituents
        FAIL,   // fail the molecule (this is equivalent to having the requireLabels option is set 'true').
    };

    struct RGoupDecompositionOptions {
        bool LabelledCores;         //allow cores to have labelled attachment points 
        bool RequireLabels;         //only allow substitutions at labelled points
        bool Symmetrize;            //attempt to symmetrize the R group assignments (experimental)
        bool KeepNonmatchingMols;   //retain molecules that do not match any of the known cores, or that have disallowed substituents, in the output.
        bool RejectDoubleAttachments;    //do not reject molecules where sidechains end up with two attachment points
        RGoupDecompositionSubstituentHandling NonLabelledSubstituentHandling; //a more flexible way of specifying what to do with substituents 
                                         // in non-labelled positions when labelledCores is true. 
        bool Verbose;
    public:
        RGoupDecompositionOptions() : LabelledCores(false)
            , RequireLabels(false)
            , Symmetrize(false)
            , KeepNonmatchingMols(false)
            , RejectDoubleAttachments(true)
            , NonLabelledSubstituentHandling(PASS)
            , Verbose(false)
            {}
    };

/* Some problems with the current python code that should be fixed in the rewrite:
•	it should support two R groups on a single center
•	it should properly support substituents at chiral centers
•	it should support multiple cores (the first one that matches a molecule is used)
*/

  //! \brief Returns ...
  //!
  /*!
      \param mols       the set of input moplecules
      \param options    options
      \param results    returned results
  */
    void RGroupDecomposite(const std::vector<ROMOL_SPTR> &mols, const std::vector<ROMOL_SPTR> &cores, const RGoupDecompositionOptions &options, std::vector<ROMOL_SPTR> &results);
 
}

#endif




