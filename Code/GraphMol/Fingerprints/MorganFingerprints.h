//
//
//  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
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
//  Created by Greg Landrum, July 2008
//
//

/*! \file MorganFingerprints.h

*/
#ifndef __RD_MORGANFPS_H__
#define __RD_MORGANFPS_H__

#include <vector>
#include <DataStructs/SparseIntVect.h>
#include <boost/cstdint.hpp>

namespace RDKit {
  class ROMol;
  namespace MorganFingerprints {
    const std::string morganFingerprintVersion="0.1.0";
    
    //! returns the Morgan fingerprint for a molecule
    /*!  
      These fingerprints are similar to the well-known ECFP or
      FCFP fingerprints, depending on which invariants are used.
        
      The algorithm used is described in the paper
      D. Rogers, R.D. Brown, M. Hahn J. Biomol. Screen. 10:682-6 (2005)
      and in more detail in an unpublished technical report:
      http://www.ics.uci.edu/~welling/teaching/ICS274Bspring06/David%20Rogers%20-%20ECFP%20Manuscript.doc

      \param mol:    the molecule to be fingerprinted
      \param radius: the number of iterations to grow the fingerprint
      \param invariants : optional pointer to a set of atom invariants to
            be used. By default ECFP-type invariants are used 
            (calculated by getConnectivityInvariants())

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    SparseIntVect<boost::uint32_t> *
      getFingerprint(const ROMol &mol,
                     unsigned int radius,
                     std::vector<boost::uint32_t> *invariants=0);
      
    //! returns the connectivity invariants for a molecule
    /*!  

      \param mol:    the molecule to be considered
      \param invariants : used to return the results
      \param includeRingMembership : if set, whether or not the atom is in
                 a ring will be used in the invariant list.
    */
    void getConnectivityInvariants(const ROMol &mol,
                                   std::vector<boost::uint32_t> &invars,
                                   bool includeRingMembership=false);

  } // end of namespace MorganFingerprints
}

#endif
