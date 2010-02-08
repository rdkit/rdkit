#include "RDKitImpl.h"
#include "RDMolImpl.h"
#include "RDMolSupplierImpl.h"
#include "nsIGenericFactory.h"
#include "nsISupportsUtils.h"

NS_GENERIC_FACTORY_CONSTRUCTOR(RDKitImpl);
NS_GENERIC_FACTORY_CONSTRUCTOR(RDMolecule);
NS_GENERIC_FACTORY_CONSTRUCTOR(RDMolSupplier);

static const nsModuleComponentInfo components[] = 
{ 
  { "RDKit Interface",  
    IRDKIT_IID,
    "@rationaldiscovery.com/RDKit/base",
    RDKitImplConstructor
  },
  { "RDKit Molecule Interface",  
    IRDMOLECULE_IID,
    "@rationaldiscovery.com/RDKit/molecule",
    RDMoleculeConstructor
  },
  { "RDKit Molecule Supplier Interface",  
    IRDMOLSUPPLIER_IID,
    "@rationaldiscovery.com/RDKit/molsupplier",
    RDMolSupplierConstructor
  },
  
}; 
 
NS_IMPL_NSGETMODULE(nsRDKitModule, components);



