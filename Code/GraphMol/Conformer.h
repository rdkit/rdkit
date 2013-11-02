//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_CONFORMER_H
#define _RD_CONFORMER_H

#include <Geometry/point.h>
#include <RDGeneral/types.h>
#include <boost/smart_ptr.hpp>

namespace RDKit {
  class ROMol;
  
  //! used to indicate errors from incorrect confomer access
  class ConformerException : public std::exception {
  public:
    //! construct with an error message
    ConformerException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    ConformerException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~ConformerException () throw () {};
  private:
    std::string _msg;
  };


  //! The class for representing 2D or 3D conformation of a molecule
  /*!
    This class contains
    - a pointer to the owing molecule
    - a vector of 3D points (positions of atoms)
  */
  class Conformer {
  public:

    friend class ROMol;

    //! Constructor
    Conformer() : df_is3D(true), d_id(0), dp_mol(NULL) {
      d_positions.clear();
    };

    //! Constructor with number of atoms specified ID specification
    Conformer(unsigned int numAtoms) : df_is3D(true), d_id(0), dp_mol(NULL) {
      if(numAtoms){
        d_positions.resize(numAtoms, RDGeom::Point3D(0.0, 0.0, 0.0));
      } else {
        d_positions.resize(0);
        d_positions.clear();
      }
    };

    //! Copy Constructor: initialize from a second conformation.
    Conformer(const Conformer &other);

    //! Destructor
    ~Conformer() {};

    //! Resize the conformer so that more atoms location can be added.
    //! Useful, for e.g., when adding hydrogens 
    void resize(unsigned int size) {
      d_positions.resize(size);
    }

    //! Reserve more space for atom position 
    void reserve(unsigned int size) {
      d_positions.reserve(size);
    }

    //! Get the molecule that oqns this conformation
    ROMol &getOwningMol() const {return *dp_mol;}

    //! Get a const reference to the vector of atom positions
    const RDGeom::POINT3D_VECT &getPositions() const;
      
    //! Get a reference to the atom positions
    RDGeom::POINT3D_VECT &getPositions(); 

    //! Get the position of the specified atom
    const RDGeom::Point3D &getAtomPos(unsigned int atomId) const; 

    //! Get the position of the specified atom
    RDGeom::Point3D &getAtomPos(unsigned int atomId); 

    //! Set the position of the specified atom
    inline void setAtomPos(unsigned int atomId, const RDGeom::Point3D &position) {
      //RANGE_CHECK(0,atomId,d_positions.size()-1);
      if (atomId >= d_positions.size()) {
        d_positions.resize(atomId+1, RDGeom::Point3D(0.0, 0.0, 0.0));
      }
      d_positions[atomId] = position;
    }

    //! get the ID of this conformer
    inline unsigned int getId() const {return d_id;}
    
    //! set the ID of this conformer
    inline void setId(unsigned int id) {
      d_id = id;
    }

    //! Get the number of atoms
    inline unsigned int getNumAtoms() const {
      return d_positions.size();
    }

    inline bool is3D() const {
      return df_is3D;
    }
    inline void set3D(bool v) {
      df_is3D=v;
    }
  protected:
    //! Set owning moelcule
    void setOwningMol(ROMol *mol);

    //! Set owning moelcule
    void setOwningMol(ROMol &mol);

  private:
    bool df_is3D; // is this a 3D conformation?
    unsigned int d_id; // id is the conformation
    ROMol *dp_mol; // owning molecule
    RDGeom::POINT3D_VECT d_positions; // positions of the atoms
  };

  typedef boost::shared_ptr<Conformer>    CONFORMER_SPTR;
}

#endif
