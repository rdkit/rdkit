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

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/StereoGroup.h>

namespace RDKit::EnumerateStereoisomers::details {
// The classes are called Flipper because that's what they were
// called in the Python, but in fact they set a particular
// stereochemistry at a centre according to whether a bool is
// true or not.
struct Flipper {
 public:
  Flipper() = default;
  Flipper(const Flipper &other) = delete;
  Flipper(Flipper &&other) = delete;
  Flipper &operator=(const Flipper &other) = delete;
  Flipper &operator=(Flipper &&other) = delete;
  virtual ~Flipper() = default;

  // Sets the stereo at the centre, one way if flag is true,
  // the other way if it is false.
  virtual void flip(bool flag) = 0;
};

struct AtomFlipper : public Flipper {
 public:
  AtomFlipper() = delete;
  AtomFlipper(RWMol &mol, const Chirality::StereoInfo &si);
  AtomFlipper(const AtomFlipper &other) = delete;
  AtomFlipper(AtomFlipper &&other) = delete;
  AtomFlipper &operator=(const AtomFlipper &other) = delete;
  AtomFlipper &operator=(AtomFlipper &&other) = delete;
  ~AtomFlipper() override = default;

  void flip(bool flag) override;

  Atom *dp_atom{nullptr};
};

struct BondFlipper : public Flipper {
 public:
  BondFlipper() = delete;
  // This c'tor may leave dp_bond as a nullptr if the bond
  // is not appropriate.
  BondFlipper(RWMol &mol, const Chirality::StereoInfo &si);
  BondFlipper(const BondFlipper &other) = delete;
  BondFlipper(BondFlipper &&other) = delete;
  BondFlipper &operator=(const BondFlipper &other) = delete;
  BondFlipper &operator=(BondFlipper &&other) = delete;
  ~BondFlipper() override = default;

  void flip(bool flag) override;

  Bond *dp_bond{nullptr};
};

struct StereoGroupFlipper : public Flipper {
 public:
  StereoGroupFlipper() = delete;
  StereoGroupFlipper(const StereoGroup &sg);
  StereoGroupFlipper(const StereoGroupFlipper &other) = delete;
  StereoGroupFlipper(StereoGroupFlipper &&other) = delete;
  StereoGroupFlipper &operator=(const StereoGroupFlipper &other) = delete;
  StereoGroupFlipper &operator=(StereoGroupFlipper &&other) = delete;
  ~StereoGroupFlipper() override = default;

  void flip(bool flag) override;

  std::vector<std::pair<Atom *, Atom::ChiralType>> d_original_parities;
};

struct AtropisomerFlipper : public Flipper {
 public:
  AtropisomerFlipper() = delete;
  AtropisomerFlipper(RWMol &mol, const Chirality::StereoInfo &si);
  AtropisomerFlipper(const AtropisomerFlipper &other) = delete;
  AtropisomerFlipper(AtropisomerFlipper &&other) = delete;
  AtropisomerFlipper &operator=(const AtropisomerFlipper &other) = delete;
  AtropisomerFlipper &operator=(AtropisomerFlipper &&other) = delete;
  ~AtropisomerFlipper() override = default;

  void flip(bool flag) override;

  Bond *dp_bond{nullptr};
  std::vector<unsigned int> d_ctrlAtoms;
};

}  // namespace RDKit::EnumerateStereoisomers::details

#endif  // FLIPPERS_H
