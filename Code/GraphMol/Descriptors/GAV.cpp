#include <iostream>
#include <vector>
#include <map>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <ForceField/MMFF/Params.h>
#include <boost/thread.hpp>
#include <boost/foreach.hpp>



/*
export BOOST=/usr/local/Cellar/boost/1.61.0/include/
export RDBASEINCLUDE=/usr/local/Cellar/rdkit/2016.03.3/include/rdkit/
export RDBASELIB=/usr/local/Cellar/rdkit/2016.03.3/lib
export BOOSTINCLUDE=/usr/local/Cellar/boost/1.61.0/include/boost/
export BOOSTLIB=/usr/local/Cellar/boost/1.61.0/lib/
clang++ -g -Wall -I$RDBASEINCLUDE -I$BOOSTINCLUDE GAV_heatComb.cpp -o GAV_heatComb  -L$RDBASELIB -L$BOOSTLIB -lRDGeneral -lGraphMol -lSmilesParse  -lFileParsers -lSubstructMatch -lForceField -lForceFieldHelpers -lboost_thread-mt -lboost_system-mt
*/

using namespace RDKit;


bool isRingAromatic(const ROMol *mol, std::vector<int> BondRing){
    int nb = BondRing.size();
        for (unsigned int i = 0; i < nb; ++i) {
          if (!mol->getBondWithIdx(BondRing[i])->getIsAromatic()){
            return false;
          }
        }
      return true;
    }


std::vector<MatchVectType> aromatic_6_NC_atoms(const ROMol *mol){
    ROMol *r6=SmartsToMol("[#6,#7;a]1[#6,#7;a][#6,#7;a][#6,#7;a][#6,#7;a][#6,#7;a]1");
    std::vector<MatchVectType> atomkeep;
    unsigned int nMatches;
    nMatches=SubstructMatch(*mol,*r6,atomkeep);
    return atomkeep;
}


VECT_INT_VECT aromatic_6_NC_bond(const ROMol *mol){
    VECT_INT_VECT aromatic6bondRings;
    VECT_INT_VECT bondRings = mol->getRingInfo()->bondRings();
    int nb= bondRings.size();
    for (unsigned int i = 0; i < nb; ++i) {
        if ((bondRings[i].size()==6) and (isRingAromatic(mol,bondRings[i]))) {
            aromatic6bondRings.push_back(bondRings[i]);

        }
    }
    return aromatic6bondRings;
}

void Aromatics6ring2(RWMol &mol,std::vector<MatchVectType> atomkeep, VECT_INT_VECT bondRings) {
   // # greg version
  // iter other vector of vector double loop!
  for (std::vector<MatchVectType>::const_iterator mvIt = atomkeep.begin();
         mvIt != atomkeep.end(); ++mvIt) {
      for (MatchVectType::const_iterator mIt = mvIt->begin();
           mIt != mvIt->end(); ++mIt) {
               mol.getAtomWithIdx(mIt->second)->setIsAromatic(true);
      }
  }

  Bond *bnd;
  int nb= bondRings.size();
  for (unsigned int i = 0; i < nb; ++i) {
      std::vector<int> bondring = bondRings[i];
      int nc = bondring.size();
      for (unsigned int j = 0; j < nc; ++j) {

               bnd = mol.getBondWithIdx(bondring[j]);
               int atombeginidx=mol.getAtomWithIdx(bnd->getBeginAtomIdx())->getAtomicNum();
               int atomendidx = mol.getAtomWithIdx(bnd->getEndAtomIdx())->getAtomicNum();

               if ((atombeginidx ==6 or atombeginidx ==7 ) and (atomendidx ==6 or atomendidx ==7)) {
               bnd->setBondType(Bond::AROMATIC);
               bnd->setIsAromatic(true);
            }
      }
  }
}

// generate the Frags like Naef:
std::string iterFrags2(ROMol *mol){
  INT_VECT allAtms;

   unsigned int numatom = mol->getNumAtoms(true);
   unsigned int nbnds = mol->getNumBonds();
   // example of bonds iterator!
   /*for (ROMol::BondIterator bi = mol->beginBonds(); bi != mol->endBonds();
         ++bi) {
      (*bi)->setIsAromatic(false);
    }
    */

        std::string frag="";

    int na = 0;
    for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms(); ++ai) {
           frag="";

        // iter other the neighbors
        int center = (*ai)->getIdx();

        if (((*ai)->getAtomicNum() == 5  or (*ai)->getAtomicNum() == 6  or
             (*ai)->getAtomicNum() == 7  or (*ai)->getAtomicNum() == 8  or
             (*ai)->getAtomicNum() == 14 or (*ai)->getAtomicNum() == 15 or
             (*ai)->getAtomicNum() == 16 ) and numatom > 1 ){
        
          na+=1;
          frag += (*ai)->getSymbol();

          Atom::HybridizationType hyb = (*ai)->getHybridization();
          if      (hyb == RDKit::Atom::SP3)  {  frag += " sp3"; }
          else if (hyb == RDKit::Atom::SP2)  {  frag += " sp2"; }
          else if (hyb == RDKit::Atom::SP)   {  frag += " sp";  }
          else  frag +=" ?";

          

            std::vector<int> neighbours;
            ROMol::ADJ_ITER nbrIdx, endNbrs;
            int nh=0;
            int nv=0;
            boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(*ai);
            while (nbrIdx != endNbrs) {

              int numatom=mol->getAtomWithIdx(*nbrIdx)->getAtomicNum();
                if (numatom ==1) nh += 1;
                if (numatom > 1) nv+=1;

                  neighbours.push_back(*nbrIdx);
           ++nbrIdx;
          }


      std::cout << frag << ":" << nh << "," << nv <<"\n";


      }

       /*
           
            nh=0
            Bondsymbols=[]
            Atomsymbols=[]
            Pisymbols=[]
            Atomnums = []
            Aros=[]
            Vals=0
            nv=0

            centersinglebond=singlebonds(atom,center,mol)

            for nei in neighbours:
                numatom=mol.GetAtomWithIdx(nei).GetAtomicNum()
                if (numatom ==1):
                    nh += 1
                if (numatom > 1):
                    nv+=1
                    bondsymbol,val,single,aro= BT(mol.GetBondBetweenAtoms(nei, center).GetBondType())
                    # count number of single bonds with heavy atoms
                    Vals+=single
                    Bondsymbols.append(bondsymbol)
                    Aros.append(aro)
                    Atomsymbols.append(re.sub('[\[\]]', '',mol.GetAtomWithIdx(nei).GetSmarts()).title())
                    if atom.GetAtomicNum() in [7,8,16]:
                        c2pi= C2(mol,center,nei)
                    else:
                        c2pi=0
                    Pisymbols.append(PI(centersinglebond,c2pi,atom.GetAtomicNum())) # Utils.NumPiElectrons(mol.GetAtomWithIdx(nei)) ; mol.GetBondBetweenAtoms(nei,center).GetIsConjugated()
                    Atomnums.append(10*numatom+val) # determine the order of the atom list using numatom + bondtypte_number
            AtomCode = []
            # determine number of H
            atomTotalH= atom.GetNumExplicitHs() + atom.GetNumImplicitHs() + nh
            # determine pi case ...  not completly correct!
            pisymbol = sum(Pisymbols)
            px=''
            if pisymbol>1:
                px='('+str(pisymbol)+'Pi)'
            elif pisymbol==1:
                px='(Pi)'
            for bondsymbol,atomsymbol in zip(Bondsymbols,Atomsymbols):
                AtomCode.append(bondsymbol+atomsymbol)
            AtomCode= [x for (y,x) in sorted(zip(Atomnums,AtomCode))]
            if len(AtomCode)+atomTotalH>=2 and nv>0: # and na>1 not good for enumerate! only look at two atoms neighbours
                AC={}
                for i in AtomCode: # create a dictionary (key,val)
                    if not AC.has_key(i): AC[i]=1
                    else: AC[i]+=1
                Code=" "
                if atomTotalH>1:
                    Code+="H"+str(atomTotalH)
                elif atomTotalH==1:
                    Code+="H"
                for key in AC:
                    if AC[key]>1:
                        Code+=key+str(AC[key])
                    elif AC[key]==1:
                        Code+=key
                if any(Aros):
                    r=" aromatic"
                mycode = re.sub('[\[\]]', '', atom.GetSmarts().title())+r+Code+px
                frag.append(mycode)
*/
      }

        if (na<2) { return ""; }
    return frag;
}





int main (int argc, char **argv)
{
  bool showHelp = (argc < 2);
  std::string smi;
  if (!showHelp) {
    smi = argv[1];
    showHelp = ((smi == "-h") || (smi == "-H") || (smi == "--help"));
  }
  if (showHelp) {
    std::cerr
      << "./GAV_heatComb <molfile>" << std::endl;
    exit(0);
  }

  // keep the H from the sdf => false is used!
  SDMolSupplier supplier(smi,false,false);

  int i = 0;

while (!supplier.atEnd())
   {
   ROMol *mol = supplier.next();
   i++;

  // get the 6 rings aromatic bonds (onyl N or C)
  VECT_INT_VECT aromaticbondRings = aromatic_6_NC_bond(mol);

  // get the 6 rings aromatic atoms (only N or C)
  std::vector<MatchVectType> matchVect = aromatic_6_NC_atoms(mol);

  // cast ROMol to RWMol
  RWMol &wmol = static_cast<RWMol &>(*mol);

  // kekulize
  MolOps::Kekulize(wmol, true);

  // rearomatize bonds & atoms
  Aromatics6ring2(wmol, matchVect, aromaticbondRings);

  // recast to ROMol
  ROMol *remol = static_cast<ROMol *>(&wmol);

  // convert to smile
  std::string smi1 = MolToSmiles(*remol, true);

  std::cout << "before sanitizeMol: " << smi1 << "\n";


  //MolOps::sanitizeMol(wmol);

  // optimize
  //MMFF::MMFFOptimizeMolecule(wmol);


  iterFrags2(mol);

  // recast to ROMol
  ROMol *remol2 = static_cast<ROMol *>(&wmol);

  // convert to smile
  std::string smi2 = MolToSmiles(*remol2, true);

  // print result
  std::cout << "result: " << smi2 << "\n";




  //std::cout << MolToMolBlock(*remol, true, -1) << std::endl;

   }
  return 0;
}
