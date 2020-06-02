//
//  Copyright (C) 2009 Greg Landrum
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_RXNPICKLE_H_2JUNE2009
#define RD_RXNPICKLE_H_2JUNE2009

#include <GraphMol/MolPickler.h>
// Std stuff
#include <iostream>
#include <string>
#include <exception>
#ifdef WIN32
#include <ios>
#endif

namespace RDKit {
class ChemicalReaction;

//! used to indicate exceptions whilst pickling (serializing) reactions
class RDKIT_CHEMREACTIONS_EXPORT ReactionPicklerException
    : public std::exception {
 public:
  ReactionPicklerException(const char *msg) : _msg(msg){};
  ReactionPicklerException(const std::string msg) : _msg(msg){};
  const char *what() const noexcept override { return _msg.c_str(); };
  ~ReactionPicklerException() noexcept {};

 private:
  std::string _msg;
};

//! handles pickling (serializing) reactions
class RDKIT_CHEMREACTIONS_EXPORT ReactionPickler {
 public:
  static const std::int32_t versionMajor;  //!< mark the pickle version
  static const std::int32_t versionMinor;  //!< mark the pickle version
  static const std::int32_t versionPatch;  //!< mark the pickle version
  static const std::int32_t endianId;  //! mark the endian-ness of the pickle

  //! the pickle format is tagged using these tags:
  //! NOTE: if you add to this list, be sure to put new entries AT THE BOTTOM,
  //! otherwise you will break old pickles.
  typedef enum {
    VERSION = 10000,
    BEGINREACTANTS,
    ENDREACTANTS,
    BEGINPRODUCTS,
    ENDPRODUCTS,
    BEGINAGENTS,
    ENDAGENTS,
    ENDREACTION,
    BEGINPROPS,
    ENDPROPS
  } Tags;

  //! pickles a reaction and sends the results to stream \c ss
  static void pickleReaction(const ChemicalReaction *rxn, std::ostream &ss,
                             unsigned int propertyFlags);
  static void pickleReaction(const ChemicalReaction *rxn, std::ostream &ss);
  static void pickleReaction(const ChemicalReaction &rxn, std::ostream &ss) {
    ReactionPickler::pickleReaction(&rxn, ss);
  };
  static void pickleReaction(const ChemicalReaction &rxn, std::ostream &ss,
                             unsigned int propertyFlags) {
    ReactionPickler::pickleReaction(&rxn, ss, propertyFlags);
  };
  //! pickles a reaction and adds the results to string \c res
  static void pickleReaction(const ChemicalReaction *rxn, std::string &res,
                             unsigned int propertyFlags);
  static void pickleReaction(const ChemicalReaction *rxn, std::string &res);
  static void pickleReaction(const ChemicalReaction &rxn, std::string &res) {
    ReactionPickler::pickleReaction(&rxn, res);
  };
  static void pickleReaction(const ChemicalReaction &rxn, std::string &res,
                             unsigned int propertyFlags) {
    ReactionPickler::pickleReaction(&rxn, res, propertyFlags);
  };

  //! constructs a reaction from a pickle stored in a
  //! string
  static void reactionFromPickle(const std::string &pickle,
                                 ChemicalReaction *rxn);
  static void reactionFromPickle(const std::string &pickle,
                                 ChemicalReaction &rxn) {
    ReactionPickler::reactionFromPickle(pickle, &rxn);
  };

  //! constructs a reaction from a pickle stored in a
  //! stream
  static void reactionFromPickle(std::istream &ss, ChemicalReaction *rxn);
  static void reactionFromPickle(std::istream &ss, ChemicalReaction &rxn) {
    ReactionPickler::reactionFromPickle(ss, &rxn);
  };

 private:
  //! do the actual work of pickling a reaction
  static void _pickle(const ChemicalReaction *rxn, std::ostream &ss,
                      unsigned int propertyFlags);

  //! do the actual work of de-pickling a reaction
  static void _depickle(std::istream &ss, ChemicalReaction *rxn, int version);

  //! pickle standard properties
  static void _pickleProperties(std::ostream &ss, const RDProps &props,
                                unsigned int pickleFlags);
  //! unpickle standard properties
  static void _unpickleProperties(std::istream &ss, RDProps &props);
};
};  // namespace RDKit

#endif
