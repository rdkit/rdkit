#include "RDMolImpl.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/Fileparsers.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <string>
#include <vector>
#include <nsMemory.h>

NS_IMPL_ISUPPORTS1(RDMolecule, IRDMolecule)

RDMolecule::~RDMolecule()
{
  if(this->dp_mol){
    delete dp_mol;
    dp_mol=0;
  } 
    
}

/* double GetMW (); */
NS_IMETHODIMP RDMolecule::GetMW(double *_retval)
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  RDKit::PeriodicTable *pt=RDKit::PeriodicTable::getTable();
  *_retval=0.0;
  for(RDKit::ROMol::AtomIterator atIt=dp_mol->beginAtoms();
      atIt!=dp_mol->endAtoms();atIt++){
    *_retval += (*atIt)->getMass()+(*atIt)->getNumImplicitHs()*pt->getAtomicWeight(1);
  }
  return NS_OK;
}

/* string GetSmiles (); */
NS_IMETHODIMP RDMolecule::GetSmiles(char **_retval)
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  std::string smi=RDKit::MolToSmiles(*dp_mol);
  *_retval=(char *)nsMemory::Clone(smi.c_str(),sizeof(char)*(smi.size()+1));
  return *_retval ? NS_OK : NS_ERROR_OUT_OF_MEMORY;
}

/* unsigned long GetSmartsMatchCount (in string smarts); */
NS_IMETHODIMP RDMolecule::GetSmartsMatchCount(const char *smarts, PRUint32 *_retval)
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  RDKit::ROMol *patt=RDKit::SmartsToMol(std::string(smarts));
  if(!patt) return NS_ERROR_FAILURE;

  std::vector<RDKit::MatchVectType> matches;
  int res=RDKit::SubstructMatch(dp_mol,patt,matches);
  *_retval = res;
  return NS_OK;
}

/* double LogP (); */
NS_IMETHODIMP RDMolecule::LogP(double *_retval)
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  double logp;
  if(dp_mol->hasProp("_CrippenLogP")) {
    dp_mol->getProp("_CrippenLogP",logp);
  } else {
    double mr;
    RDKit::Descriptors::CalcCrippenDescriptors(dp_mol,logp,mr);
    dp_mol->setProp("_CrippenLogP",logp,true);
    dp_mol->setProp("_CrippenMR",mr,true);
  }
  *_retval=logp;
  return NS_OK;
}

/* double MR (); */
NS_IMETHODIMP RDMolecule::MR(double *_retval)
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  double mr;
  if(dp_mol->hasProp("_CrippenMR")) {
    dp_mol->getProp("_CrippenMR",mr);
  } else {
    double logp;
    RDKit::Descriptors::CalcCrippenDescriptors(dp_mol,logp,mr);
    dp_mol->setProp("_CrippenLogP",logp,true);
    dp_mol->setProp("_CrippenMR",mr,true);
  }
  *_retval=mr;
  return NS_OK;
}

/* string GetMolBlock (); */
NS_IMETHODIMP RDMolecule::GetMolBlock(char **_retval)
{

  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  std::string molB=RDKit::MolToMolBlock(dp_mol);
  *_retval=(char *)nsMemory::Clone(molB.c_str(),sizeof(char)*(molB.size()+1));
  return *_retval ? NS_OK : NS_ERROR_OUT_OF_MEMORY;
}

/* void Generate3DCoords (); */
NS_IMETHODIMP RDMolecule::Generate3DCoords()
{
  if(!dp_mol) return NS_ERROR_NOT_INITIALIZED;
  bool embedded=RDKit::DGeomHelpers::EmbedMolecule(*dp_mol);
  if(!embedded) return NS_ERROR_FAILURE;
  
  ForceFields::ForceField *ff=RDKit::UFF::constructForceField(dp_mol);
  if(!ff) return NS_ERROR_FAILURE;
  ff->initialize();
  int needsMore=ff->minimize();
  delete ff;
  if(needsMore) return NS_ERROR_FAILURE;
  return NS_OK;
}

