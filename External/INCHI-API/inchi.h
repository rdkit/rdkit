//
//  Copyright (c) 2011-2014, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written
//       permission.
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

#include <RDGeneral/export.h>
#ifndef RDKIT_INCHI_30JUNE2011
#define RDKIT_INCHI_30JUNE2011
#include <GraphMol/RDKitBase.h>
#include <string>
namespace RDKit {
struct RDKIT_RDINCHILIB_EXPORT ExtraInchiReturnValues {
  int returnCode;
  std::string messagePtr;
  std::string logPtr;
  std::string auxInfoPtr;  // not used for InchiToMol
};
/*! Get a RWMol molecule instance for a InChI string.
 * \param inchi The input InChI string, which can be standard or not
 * \param rv An ExtraInchiReturnValues struct instance that is used to receive
 * extra return values such as error messages from InChI API.
 * \param sanitize Whether to sanitize the generated molecule before returning
 * it
 * \param removeHs Whether to remove hydrogens from the generated molecule
 * before returning it.
 */
RDKIT_RDINCHILIB_EXPORT RWMol* InchiToMol(const std::string& inchi,
                                          ExtraInchiReturnValues& rv,
                                          bool sanitize = true,
                                          bool removeHs = true);
/*! Get the InChI string for a given molecule
 * \param mol The input molecule
 * \param rv An ExtraInchiReturnValues struct instance that is used to receive
 * extra return values such as InChI Auxiliary Information and error messages
 * from InChI API.
 * \param options An null-terminated character string of space-deliminated
 * InChI options that is passed to InChI API as is (except that / is naively
 * converted to - to non-Windows platforms and - is converted to / on Windows)
 * Available options are explained in the InChI technical FAQ:
 * http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14
 * and the User Guide:
 * http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf
 */
RDKIT_RDINCHILIB_EXPORT std::string MolToInchi(const ROMol& mol,
                                               ExtraInchiReturnValues& rv,
                                               const char* options = NULL);
/*! Get the InChI string for a given mol block
 * \param mol The input mol block
 * \param rv An ExtraInchiReturnValues struct instance that is used to receive
 * extra return values such as InChI Auxiliary Information and error messages
 * from InChI API.
 * \param options An null-terminated character string of space-deliminated
 * InChI options that is passed to InChI API as is (except that / is naively
 * converted to - to non-Windows platforms and - is converted to / on Windows)
 * Available options are explained in the InChI technical FAQ:
 * http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14
 * and the User Guide:
 * http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf
 */
RDKIT_RDINCHILIB_EXPORT std::string MolBlockToInchi(const std::string& mol,
                                                    ExtraInchiReturnValues& rv,
                                                    const char* options = NULL);
/*! Get the InChI Key for an input InChI string
 * \param inchi The input InChI string, which can be standard or not.
 */
RDKIT_RDINCHILIB_EXPORT std::string InchiToInchiKey(const std::string& inchi);

/*! Get the InChI key for a given molecule directly
 * \param mol The input molecule
 * \param options An null-terminated character string of space-deliminated
 * InChI options that is passed to InChI API as is (except that / is naively
 * converted to - to non-Windows platforms and - is converted to / on Windows)
 * Available options are explained in the InChI technical FAQ:
 * http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14
 * and the User Guide:
 * http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf
 */
inline std::string MolToInchiKey(const ROMol& mol, const char* options = NULL) {
  ExtraInchiReturnValues rv;
  return InchiToInchiKey(MolToInchi(mol, rv, options));
};

}  // namespace RDKit
#endif
