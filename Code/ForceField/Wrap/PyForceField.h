// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <algorithm>
#include <Geometry/point.h>

namespace ForceFields {
  class PyForceField {
  public:
    PyForceField( ForceField *f) : field(f) {};

    ~PyForceField() {
      //std::cerr << " *** destroy PyForce field " << std::endl;
      field.reset();
      //std::cerr << " ***       reset DONE" << std::endl;
      extraPoints.clear();
      //std::cerr << " *** destroy PyForce field DONE" << std::endl;
    }
    
    int addExtraPoint(double x,double y,double z,bool fixed=true){
      RDGeom::Point3D *pt=new RDGeom::Point3D(x,y,z);
      PRECONDITION(this->field,"no force field");
      this->extraPoints.push_back(boost::shared_ptr<RDGeom::Point3D>(pt));
      unsigned int ptIdx = this->extraPoints.size()-1;
      RDGeom::Point3D *ptr=this->extraPoints[ptIdx].get();
      this->field->positions().push_back(ptr);
      int idx=this->field->positions().size();
      if(fixed){
	this->field->fixedPoints().push_back(idx-1);
      }
      return idx;
    }

    double calcEnergy() {
      PRECONDITION(this->field,"no force field");
      return this->field->calcEnergy();
    }
    
    int minimize(int maxIts,double forceTol,double energyTol){
      PRECONDITION(this->field,"no force field");
      return this->field->minimize(maxIts,forceTol,energyTol);
    }

    void initialize() {
      PRECONDITION(this->field,"no force field");
      this->field->initialize();
    }
    
    //private:
    std::vector< boost::shared_ptr<RDGeom::Point3D> > extraPoints;
    boost::shared_ptr<ForceField> field;
  };

  class PyMMFFMolProperties {
  public:
    PyMMFFMolProperties(RDKit::MMFF::MMFFMolProperties *mp) : mmffMolProperties(mp) {};
    ~PyMMFFMolProperties() {};
      
    unsigned int getMMFFAtomType(unsigned int idx) {
      return (unsigned int)(mmffMolProperties->getMMFFAtomType(idx));
    };
    double getMMFFFormalCharge(unsigned int idx) {
      return mmffMolProperties->getMMFFFormalCharge(idx);
    };
    double getMMFFPartialCharge(unsigned int idx) {
      return mmffMolProperties->getMMFFPartialCharge(idx);
    };
    void setMMFFDielectricModel(boost::uint8_t dielModel)
    {
      mmffMolProperties->setMMFFDielectricModel(dielModel);
    };
    void setMMFFDielectricConstant(double dielConst)
    {
      mmffMolProperties->setMMFFDielectricConstant(dielConst);
    };
    void setMMFFBondTerm(bool state)
    {
      mmffMolProperties->setMMFFBondTerm(state);
    };
    void setMMFFAngleTerm(const bool state)
    {
      mmffMolProperties->setMMFFAngleTerm(state);
    };
    void setMMFFStretchBendTerm(const bool state)
    {
      mmffMolProperties->setMMFFStretchBendTerm(state);
    };
    void setMMFFOopTerm(const bool state)
    {
      mmffMolProperties->setMMFFOopTerm(state);
    };
    void setMMFFTorsionTerm(const bool state)
    {
      mmffMolProperties->setMMFFTorsionTerm(state);
    };
    void setMMFFVdWTerm(const bool state)
    {
      mmffMolProperties->setMMFFVdWTerm(state);
    };
    void setMMFFEleTerm(const bool state)
    {
      mmffMolProperties->setMMFFEleTerm(state);
    };
    void setMMFFVariant(std::string mmffVariant)
    {
      mmffMolProperties->setMMFFVariant(mmffVariant);
    };
    void setMMFFVerbosity(unsigned int verbosity)
    {
      mmffMolProperties->setMMFFVerbosity(verbosity);
    };
    boost::shared_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
  };
}
