//
//  Copyright (C) 2020 Guillaume GODIN
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Adding ATOM FEATURES descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <iostream>
#include "Augmentation.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GraphMol.h>
#include <boost/foreach.hpp>
#include <cmath>
#include <vector>


namespace RDKit {

namespace Descriptors {

namespace {



void AugmentationVector(const ROMol* mol, std::vector <std::string> &res, unsigned int naug, bool addcanonical) {

    PRECONDITION(mol,"bad mol");

    if (addcanonical) {
    	res.push_back( MolToSmiles( *mol ) );
    }

    std::vector<std::string> nonuniqueres(naug);
    std::string smile;
    for (unsigned int i = 0 ; i < naug; i++) {
        smile = MolToSmiles( *mol, true, false, -1, false, false, false, true ) ;
        nonuniqueres[i] = smile;
    }

    // make a unique list
    std::sort(nonuniqueres.begin(), nonuniqueres.end());
    std::unique_copy(nonuniqueres.begin(), nonuniqueres.end(), std::back_inserter(res));
    //std::set_difference(nonuniqueres.begin(), nonuniqueres.end(), res.begin(), res.end(), std::ostream_iterator<std::string>(std::cout, " "));


}
}  // end of anonymous namespace

// entry point
void AugmentationVect(const ROMol& mol, std::vector <std::string>& res, unsigned int naug, bool addcanonical) {

  AugmentationVector( &mol, res, naug, addcanonical);
}

}  // namespace Descriptors
}  // namespace RDKit