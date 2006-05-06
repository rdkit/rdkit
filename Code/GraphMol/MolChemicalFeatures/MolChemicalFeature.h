//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
    MolChemicalFeature(const ROMol *mol,const MolChemicalFeatureFactory *factory,
                       const MolChemicalFeatureDef *fdef) :
      dp_mol(mol), dp_factory(factory), dp_def(fdef), d_activeConf(-1) {};

    ~MolChemicalFeature() {}

    //! Functions we definitely need to implement from the base class
    const std::string &getFamily() const;
    const std::string &getType() const;
    RDGeom::Point3D getPos() const;

    //! Functions that are more specific to the MolChemicalFeature
    RDGeom::Point3D getPos(int confId) const;
    const MolChemicalFeatureFactory *getFactory() const { return dp_factory; };
    const ROMol *getMol() const { return dp_mol; };
    const MolChemicalFeatureDef *getFeatDef() const { return dp_def; };
    
    inline unsigned int getNumAtoms() const {
      return d_atoms.size();
    }
    AtomPtrContainer::const_iterator beginAtoms() const { return d_atoms.begin(); };
    AtomPtrContainer::const_iterator endAtoms() const { return d_atoms.end(); };

    void setActiveConformer(int confId);

    int getActiveConformer() const { return d_activeConf;};

    void clearCache() { d_locs.clear(); };
    
    const AtomPtrContainer &getAtoms() const {
      return d_atoms;
    }

  private:
    typedef std::map<int,RDGeom::Point3D> PointCacheType;
    const ROMol *dp_mol;
    const MolChemicalFeatureFactory *dp_factory;
    const MolChemicalFeatureDef *dp_def;
    //std::string d_type;
    //std::string d_family;

    int d_activeConf;
    AtomPtrContainer d_atoms;
    mutable PointCacheType d_locs;
  };
  
}

#endif
