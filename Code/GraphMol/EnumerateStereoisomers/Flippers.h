//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef FLIPPERS_H
#define FLIPPERS_H

#include <RDGeneral/export.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/StereoGroup.h>

#include <limits>

namespace RDKit {

class ROMol;

namespace EnumerateStereoisomers {
namespace details {
// The classes are called Flipper because that's what they were
// called in the Python, but in fact they set a particular
// stereochemistry at a centre according to whether a bool is
// true or not.  Setters would be a more accurate name.
struct RDKIT_ENUMERATESTEREOISOMERS_EXPORT Flipper {
 public:
  Flipper() = default;
  Flipper(const Flipper &other) = delete;
  Flipper(Flipper &&other) = delete;
  Flipper &operator=(const Flipper &other) = delete;
  Flipper &operator=(Flipper &&other) = delete;
  virtual ~Flipper() = default;

  // Sets the stereo at the centre, one way if flag is true,
  // the other way if it is false.
  virtual void flip(ROMol &mol, bool flag) = 0;
};

class RDKIT_ENUMERATESTEREOISOMERS_EXPORT AtomFlipper : public Flipper {
 public:
  AtomFlipper() = delete;
  explicit AtomFlipper(const Chirality::StereoInfo &si);
  AtomFlipper(const AtomFlipper &other) = delete;
  AtomFlipper(AtomFlipper &&other) = delete;
  AtomFlipper &operator=(const AtomFlipper &other) = delete;
  AtomFlipper &operator=(AtomFlipper &&other) = delete;
  ~AtomFlipper() override = default;

  void flip(ROMol &mol, bool flag) override;

 private:
  unsigned int d_atomIdx{std::numeric_limits<unsigned int>::max()};
};

class RDKIT_ENUMERATESTEREOISOMERS_EXPORT BondFlipper : public Flipper {
 public:
  BondFlipper() = delete;
  // This c'tor may mark the flipper as invalid if the bond
  // is not appropriate.
  BondFlipper(ROMol &mol, const Chirality::StereoInfo &si);
  BondFlipper(const BondFlipper &other) = delete;
  BondFlipper(BondFlipper &&other) = delete;
  BondFlipper &operator=(const BondFlipper &other) = delete;
  BondFlipper &operator=(BondFlipper &&other) = delete;
  ~BondFlipper() override = default;

  bool isValid() const;
  void flip(ROMol &mol, bool flag) override;

 private:
  unsigned int d_bondIdx{std::numeric_limits<unsigned int>::max()};
  bool d_valid{true};
};

class RDKIT_ENUMERATESTEREOISOMERS_EXPORT StereoGroupFlipper : public Flipper {
 public:
  StereoGroupFlipper() = delete;
  StereoGroupFlipper(const StereoGroup &sg);
  StereoGroupFlipper(const StereoGroupFlipper &other) = delete;
  StereoGroupFlipper(StereoGroupFlipper &&other) = delete;
  StereoGroupFlipper &operator=(const StereoGroupFlipper &other) = delete;
  StereoGroupFlipper &operator=(StereoGroupFlipper &&other) = delete;
  ~StereoGroupFlipper() override = default;

  void flip(ROMol &mol, bool flag) override;

 private:
  std::vector<std::pair<unsigned int, Atom::ChiralType>> d_original_parities;
};

class RDKIT_ENUMERATESTEREOISOMERS_EXPORT AtropisomerFlipper : public Flipper {
 public:
  AtropisomerFlipper() = delete;
  explicit AtropisomerFlipper(const Chirality::StereoInfo &si);
  AtropisomerFlipper(const AtropisomerFlipper &other) = delete;
  AtropisomerFlipper(AtropisomerFlipper &&other) = delete;
  AtropisomerFlipper &operator=(const AtropisomerFlipper &other) = delete;
  AtropisomerFlipper &operator=(AtropisomerFlipper &&other) = delete;
  ~AtropisomerFlipper() override = default;

  void flip(ROMol &mol, bool flag) override;

 private:
  unsigned int d_bondIdx{std::numeric_limits<unsigned int>::max()};
  std::vector<unsigned int> d_ctrlAtoms;
};

}  // namespace details
}  // namespace EnumerateStereoisomers
}  // namespace RDKit

#endif  // FLIPPERS_H
