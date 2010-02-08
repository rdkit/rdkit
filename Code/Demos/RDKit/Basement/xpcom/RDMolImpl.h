#ifndef _RDMOLIMPL_H_
#define _RDMOLIMPL_H_
#include "IRDKit.h"

namespace RDKit {
  class ROMol;
}

class RDMolecule : public IRDMolecule
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDMOLECULE

  RDMolecule() : dp_mol(0) {};
  RDMolecule(RDKit::ROMol *mol) : dp_mol(mol) {};
  RDKit::ROMol *dp_mol;
private:
  ~RDMolecule();

protected:
  /* additional members */
};

#endif
