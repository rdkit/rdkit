#include "RDMolSupplierImpl.h"
#include "RDMolImpl.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <string>
#include <nsMemory.h>

NS_IMPL_ISUPPORTS1(RDMolSupplier, IRDMolSupplier)

RDMolSupplier::~RDMolSupplier()
{
  if(this->dp_suppl){
    delete dp_suppl;
    dp_suppl=0;
  } 
    
}

/* boolean atEnd (); */
NS_IMETHODIMP RDMolSupplier::AtEnd(PRBool *_retval)
{
  if(!dp_suppl) return NS_ERROR_NOT_INITIALIZED;
  *_retval=dp_suppl->atEnd();
  return NS_OK;
}

/* IRDMolecule next (); */
NS_IMETHODIMP RDMolSupplier::Next(IRDMolecule **_retval)
{
  if(!dp_suppl) return NS_ERROR_NOT_INITIALIZED;
  if(dp_suppl->atEnd()) return NS_BASE_STREAM_CLOSED;
  RDKit::ROMol *roMol=dp_suppl->next();
  if(!roMol) return NS_ERROR_UNEXPECTED;

  RDMolecule *mol = new RDMolecule(roMol);
  if(!mol) return NS_ERROR_OUT_OF_MEMORY;
  *_retval = static_cast<IRDMolecule *>(mol);

  // FIX: does this leak?
  NS_ADDREF(*_retval);

  return NS_OK;
}
