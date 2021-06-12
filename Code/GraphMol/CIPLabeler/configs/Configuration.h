//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <memory>
#include <vector>

#include "../Descriptor.h"
#include "../Digraph.h"
#include "../CIPMol.h"

namespace RDKit {

class Atom;
class Bond;

namespace CIPLabeler {

class Rules;

class Configuration {
 public:
  template <typename T>
  static int parity4(const std::vector<T> &trg, const std::vector<T> &ref) {
    if (ref.size() != 4 || trg.size() != ref.size()) {
      throw std::runtime_error("Parity vectors must have size 4.");
    }

    if (ref[0] == trg[0]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> a,b,c,d
        if (ref[2] == trg[2] && ref[3] == trg[3]) {
          return 2;
        }
        // a,b,c,d -> a,b,d,c
        if (ref[2] == trg[3] && ref[3] == trg[2]) {
          return 1;
        }
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> a,c,b,d
        if (ref[2] == trg[1] && ref[3] == trg[3]) {
          return 1;
        }
        // a,b,c,d -> a,c,d,b
        if (ref[2] == trg[3] && ref[3] == trg[1]) {
          return 2;
        }
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> a,d,c,b
        if (ref[2] == trg[2] && ref[3] == trg[1]) {
          return 1;
        }
        // a,b,c,d -> a,d,b,c
        if (ref[2] == trg[1] && ref[3] == trg[2]) {
          return 2;
        }
      }
    } else if (ref[0] == trg[1]) {
      if (ref[1] == trg[0]) {
        // a,b,c,d -> b,a,c,d
        if (ref[2] == trg[2] && ref[3] == trg[3]) {
          return 1;
        }
        // a,b,c,d -> b,a,d,c
        if (ref[2] == trg[3] && ref[3] == trg[2]) {
          return 2;
        }
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> b,c,a,d
        if (ref[2] == trg[0] && ref[3] == trg[3]) {
          return 2;
        }
        // a,b,c,d -> b,c,d,a
        if (ref[2] == trg[3] && ref[3] == trg[0]) {
          return 1;
        }
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> b,d,c,a
        if (ref[2] == trg[2] && ref[3] == trg[0]) {
          return 2;
        }
        // a,b,c,d -> b,d,a,c
        if (ref[2] == trg[0] && ref[3] == trg[2]) {
          return 1;
        }
      }
    } else if (ref[0] == trg[2]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> c,b,a,d
        if (ref[2] == trg[0] && ref[3] == trg[3]) {
          return 1;
        }
        // a,b,c,d -> c,b,d,a
        if (ref[2] == trg[3] && ref[3] == trg[0]) {
          return 2;
        }
      } else if (ref[1] == trg[0]) {
        // a,b,c,d -> c,a,b,d
        if (ref[2] == trg[1] && ref[3] == trg[3]) {
          return 2;
        }
        // a,b,c,d -> c,a,d,b
        if (ref[2] == trg[3] && ref[3] == trg[1]) {
          return 1;
        }
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> c,d,a,b
        if (ref[2] == trg[0] && ref[3] == trg[1]) {
          return 2;
        }
        // a,b,c,d -> c,d,b,a
        if (ref[2] == trg[1] && ref[3] == trg[0]) {
          return 1;
        }
      }
    } else if (ref[0] == trg[3]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> d,b,c,a
        if (ref[2] == trg[2] && ref[3] == trg[0]) {
          return 1;
        }
        // a,b,c,d -> d,b,a,c
        if (ref[2] == trg[0] && ref[3] == trg[2]) {
          return 2;
        }
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> d,c,b,a
        if (ref[2] == trg[1] && ref[3] == trg[0]) {
          return 2;
        }
        // a,b,c,d -> d,c,a,b
        if (ref[2] == trg[0] && ref[3] == trg[1]) {
          return 1;
        }
      } else if (ref[1] == trg[0]) {
        // a,b,c,d -> d,a,c,b
        if (ref[2] == trg[2] && ref[3] == trg[1]) {
          return 2;
        }
        // a,b,c,d -> d,a,b,c
        if (ref[2] == trg[1] && ref[3] == trg[2]) {
          return 1;
        }
      }
    }

    // We should never hit this, but the compiler still complains
    // about a missing return statement.
    return 0;
  }

  Configuration() = delete;

  Configuration(const CIPMol &mol, Atom *focus);

  Configuration(const CIPMol &mol, std::vector<Atom *> &&foci);

  virtual ~Configuration();

  Atom *getFocus() const;

  const std::vector<Atom *> &getFoci() const;

  const std::vector<Atom *> &getCarriers() const;

  Digraph &getDigraph();

  virtual Descriptor label(Node *node, Digraph &digraph, const Rules &comp);

  virtual Descriptor label(const Rules &comp) = 0;

  virtual void setPrimaryLabel(Descriptor desc) = 0;

 protected:
  Edge *findInternalEdge(const std::vector<Edge *> &edges, Atom *f1, Atom *f2);

  bool isInternalEdge(const Edge *edge, Atom *f1, Atom *f2);

  void removeInternalEdges(std::vector<Edge *> &edges, Atom *f1, Atom *f2);

  void setCarriers(std::vector<Atom *> &&carriers);

 private:
  /**
   * Foci are the atoms on which the configuration is based,
   * and which will carry the label. E.g., the chiral atom in
   * a tetrahedral chirality, or the bond ends in a double bond.
   */
  std::vector<Atom *> d_foci;

  /**
   * Carriers are the atoms neighboring the foci that define the
   * configuration. E.g., for a chiral atom, its four neighbors
   * define a parity; for a double bond, one neighbor on each
   * side of the bond defines it as Cis or Trans.
   */
  std::vector<Atom *> d_carriers;

  Digraph d_digraph;

};  // namespace CIPLabeler

}  // namespace CIPLabeler
}  // namespace RDKit
