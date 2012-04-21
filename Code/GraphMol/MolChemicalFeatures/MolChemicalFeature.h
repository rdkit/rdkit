//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __MOLCHEMICALFEATURE_H_11012005_1404__
#define __MOLCHEMICALFEATURE_H_11012005_1404__

#include <string>
#include <vector>
#include <map>
#include <Geometry/point.h>
#include <ChemicalFeatures/ChemicalFeature.h>

namespace RDKit {
  class ROMol;
  class Atom;
  class MolChemicalFeatureFactory;
  class MolChemicalFeatureDef;
  
  
  class MolChemicalFeature : public ChemicalFeatures::ChemicalFeature {

    friend class MolChemicalFeatureFactory;
  public:
    typedef std::vector<const Atom *> AtomPtrContainer;
    typedef AtomPtrContainer::const_iterator AtomPtrContainer_CI;

    //! Constructor
    MolChemicalFeature(const ROMol *mol,
                       const MolChemicalFeatureFactory *factory,
                       const MolChemicalFeatureDef *fdef,
                       int id=-1) :
      dp_mol(mol), dp_factory(factory), dp_def(fdef), d_id(id), d_activeConf(-1) {};

    ~MolChemicalFeature() {}

    //! \brief return the name of the feature's family
    const std::string &getFamily() const;
    //! \brief return the name of the feature's type
    const std::string &getType() const;
    //! \brief return the position of the feature (obtained from
    //! from the associated conformation
    RDGeom::Point3D getPos() const;

    //! \brief return the position of the feature (obtained from
    //! from the requested conformation from the associated molecule)
    RDGeom::Point3D getPos(int confId) const;
    //! \brief return a pointer to our feature factory
    const MolChemicalFeatureFactory *getFactory() const { return dp_factory; };
    //! \brief return a pointer to our associated molecule
    const ROMol *getMol() const { return dp_mol; };
    //! \brief return a pointer to our feature definition
    const MolChemicalFeatureDef *getFeatDef() const { return dp_def; };

    //! \brief returns the active conformer (in the associated molecule)
    const int getId() const { return d_id;};
    
    //! \brief returns the number of atoms defining the feature
    inline unsigned int getNumAtoms() const {
      return d_atoms.size();
    }

    //! \brief sets the active conformer (in the associated molecule)
    void setActiveConformer(int confId);

    //! \brief returns the active conformer (in the associated molecule)
    int getActiveConformer() const { return d_activeConf;};

    //! \brief clears out the internal position cache
    void clearCache() { d_locs.clear(); };
    
    //! \brief returns our atom container of
    const AtomPtrContainer &getAtoms() const {
      return d_atoms;
    }
    AtomPtrContainer::const_iterator beginAtoms() const { return d_atoms.begin(); };
    AtomPtrContainer::const_iterator endAtoms() const { return d_atoms.end(); };

  private:
    typedef std::map<int,RDGeom::Point3D> PointCacheType;
    
    const ROMol *dp_mol;
    const MolChemicalFeatureFactory *dp_factory;
    const MolChemicalFeatureDef *dp_def;
    int d_id;
    int d_activeConf;
    AtomPtrContainer d_atoms;
    mutable PointCacheType d_locs;
  };
  
}

#endif
