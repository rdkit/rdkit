// $Id$
//
// Copyright (C) 2011-2011 Novartis Institutes for BioMedical Research, Inc
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
namespace RDKit {
  struct ExtraInchiReturnValues {
    int returnCode;
    std::string* messagePtr;
    std::string* logPtr;
    std::string* auxInfoPtr;  // not used for InchiToMol
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
  RWMol* InchiToMol(const std::string &inchi, ExtraInchiReturnValues& rv,
                    bool sanitize=true, bool removeHs=true);
  /*! Get the InChI string for a given molecule
   * \param mol The input molecule
   * \param rv An ExtraInchiReturnValues struct instance that is used to receive
   * extra return values such as InChI Auxiliary Information and error messages
   * from InChI API.
   * \param options An null-terminated character string of space-deliminated
   * InChI options that is passed to InChI API as is (except that / is naively
   * converted to - to non-Windows platforms and - is converted to / on Windows)
   */
  std::string MolToInchi(const ROMol& mol, ExtraInchiReturnValues& rv,
                         const char *options=NULL);
  /*! Get the InChI Key for an input InChI string
   * \param inchi The input InChI string, which can be standard or not.
   */
  std::string InchiToInchiKey(const std::string &inchi);
}
