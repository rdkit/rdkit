//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_MMFFBUILDER_H_
#define _RD_MMFFBUILDER_H_

#include <vector>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/cstdint.hpp>


namespace ForceFields {
  class ForceField;
}

namespace RDKit {
  class ROMol;
  namespace MMFF {

    //! Builds and returns a MMFF force field for a molecule
    /*!
      
      \param mol              the molecule to use
      \param nonBondedThresh  the threshold to be used in adding non-bonded terms
                              to the force field. Any non-bonded contact whose current
			distance is greater than \c nonBondedThresh * the minimum value
			for that contact will not be included.
      \param confId     the optional conformer id, if this isn't provided, the molecule's
                        default confId will be used.
      \param ignoreInterfragInteractions if true, nonbonded terms will not be added between
                                         fragments

      \return the new force field. The client is responsible for free'ing this.
    */
    ForceFields::ForceField *constructForceField(ROMol &mol,
      double nonBondedThresh = 100.0, int confId = -1, bool ignoreInterfragInteractions = true);


    //! Builds and returns a MMFF force field for a molecule
    /*!
      
      \param mol        the molecule to use
      \param mmffMolProperties        pointer to a MMFFMolProperties (as obtained by a
                        call to setupMMFFForceField())
      \param nonBondedThresh  the threshold to be used in adding non-bonded terms
                        to the force field. Any non-bonded contact whose current
			distance is greater than \c nonBondedThresh * the minimum value
			for that contact will not be included.
      \param confId     the optional conformer id, if this isn't provided, the molecule's
                        default confId will be used.
      \param ignoreInterfragInteractions if true, nonbonded terms will not be added between
                                         fragments
    
      \return the new force field. The client is responsible for free'ing this.
    */
    ForceFields::ForceField *constructForceField(ROMol &mol, MMFFMolProperties *mmffMolProperties,
      double nonBondedThresh = 100.0, int confId = -1, bool ignoreInterfragInteractions = true);

    namespace Tools {
      enum {
        RELATION_1_2 = 0,
        RELATION_1_3 = 1,
        RELATION_1_4 = 2,
        RELATION_1_X = 3
      };
      // these functions are primarily exposed so they can be tested.
      void setTwoBitCell(boost::shared_array<boost::uint8_t> &res,
        unsigned int pos, boost::uint8_t value);
      boost::uint8_t getTwoBitCell
        (boost::shared_array<boost::uint8_t> &res, unsigned int pos);
      boost::shared_array<boost::uint8_t> buildNeighborMatrix(const ROMol &mol);
      void addBonds(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field);
      void addAngles(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field);
      void addStretchBend(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field);
      void addOop(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field);
      void addTorsions(const ROMol &mol,
        MMFFMolProperties *mmffMolProperties, ForceFields::ForceField *field,
        std::string torsionBondSmarts = "[!$(*#*)&!D1]~[!$(*#*)&!D1]");
      void addVdW(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
        ForceFields::ForceField *field, boost::shared_array<boost::uint8_t> neighborMatrix,
        double nonBondedThresh = 100.0, bool ignoreInterfragInteractions = true);
      void addEle(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
        ForceFields::ForceField *field, boost::shared_array<boost::uint8_t> neighborMatrix,
        double nonBondedThresh = 100.0, bool ignoreInterfragInteractions = true);
    }
  }
}


#endif
