#include "RDKitImpl.h"
#include "RDMolImpl.h"
#include "RDMolSupplierImpl.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <string>

/* Implementation file */
NS_IMPL_ISUPPORTS1(RDKitImpl, IRDKit)

RDKitImpl::RDKitImpl()
{
  /* member initializers and constructor code */
}


RDKitImpl::~RDKitImpl()
{
  /* destructor code */
}

/* unsigned long strlen (in string arg); */
NS_IMETHODIMP RDKitImpl::Strlen(const char *arg, PRUint32 *_retval)
{
  std::string text(arg);
  *_retval = text.size();
  return NS_OK;
}

/* IRDMolecule MolFromSmiles (in string smiles); */
NS_IMETHODIMP RDKitImpl::MolFromSmiles(const char *smiles, IRDMolecule **_retval)
{
  std::string smi(smiles);
  RDKit::ROMol *roMol = RDKit::SmilesToMol(smiles);
  if(!roMol) return NS_ERROR_FAILURE;

  RDMolecule *mol = new RDMolecule(roMol);
  if(!mol) return NS_ERROR_OUT_OF_MEMORY;
  *_retval = static_cast<IRDMolecule *>(mol);

  // FIX: does this leak?
  NS_ADDREF(*_retval);

  return NS_OK;
}

/* IRDMolecule MolFromMolBlock (in string molBlock); */
NS_IMETHODIMP RDKitImpl::MolFromMolBlock(const char *molBlock, IRDMolecule **_retval)
{
  RDKit::ROMol *roMol = RDKit::MolBlockToMol(std::string(molBlock));
  if(!roMol) return NS_ERROR_FAILURE;

  RDMolecule *mol = new RDMolecule(roMol);
  if(!mol) return NS_ERROR_OUT_OF_MEMORY;
  *_retval = static_cast<IRDMolecule *>(mol);

  // FIX: does this leak?
  NS_ADDREF(*_retval);

  return NS_OK;
}

  /* IRDMolSupplier SupplierFromSDFile (in string fileName); */
NS_IMETHODIMP RDKitImpl::SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval)
{
  RDKit::MolSupplier *sdSuppl = new RDKit::SDMolSupplier(std::string(fileName));
  if(!sdSuppl) return NS_ERROR_FAILURE;

  RDMolSupplier *suppl = new RDMolSupplier(sdSuppl);
  if(!suppl) return NS_ERROR_OUT_OF_MEMORY;
  *_retval = static_cast<IRDMolSupplier *>(suppl);

  // FIX: does this leak?
  NS_ADDREF(*_retval);

  return NS_OK;
}

