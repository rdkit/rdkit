#ifndef _RDMOLSUPPLIERIMPL_H_
#define _RDMOLSUPPLIERIMPL_H_
#include "IRDKit.h"

namespace RDKit {
  class ROMol;
  class MolSupplier;
}

class RDMolSupplier : public IRDMolSupplier
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDMOLSUPPLIER

  RDMolSupplier() : dp_suppl(0) {};
  RDMolSupplier(RDKit::MolSupplier *suppl) : dp_suppl(suppl) {};
  RDKit::MolSupplier *dp_suppl;
private:
  ~RDMolSupplier();

protected:
  /* additional members */
};

#endif
