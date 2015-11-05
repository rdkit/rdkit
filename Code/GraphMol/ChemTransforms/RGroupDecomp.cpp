//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>

#include "../MolOps.h"
#include "../SmilesParse/SmilesParse.h"
//#include "../SmilesParse/SmilesWrite.h"
#include "../Substruct/SubstructMatch.h"
#include "ChemTransforms.h"

#include "RGroupDecomp.h"

namespace RDKit{ 

//typedef std::vector<ROMOL_SPTR>::const_iterator MolIterator;
//TODO:                          
typedef ROMOL_SPTR sidechain_t;


static
void GetMolDegenPoints (ROMol& mol, std::vector<unsigned>& equivPoints) {
    std::vector<MatchVectType> matches;
    GetSubstructMatch(mol, mol, matches, false);
    for(unsigned mIdx=0; mIdx < matches.size(); mIdx++) {
        for(unsigned idx=0; idx < matches[mIdx].size(); idx++) {
            if(matches[mIdx][idx][0] != matches[mIdx][idx][1])
                equivPoints.push_back(matches[mIdx][idx][0]);
            equivPoints.push_back(matches[mIdx][idx][1]);
        }
    }
}

//TODO:                          
/*
def SymmetrizeSidechains(mols,core,sidechains,options=None):
  if len(mols)!=len(sidechains):
    raise ValueError('sidechains list must be as long as molecules list')
  if not core:
    return sidechains

  matches = core.GetSubstructMatches(core,uniquify=False)
  if len(matches)<=1:
    return sidechains

  degenPts = GetMolDegenPoints(core)

  # start by ordering the molecules by the number of attachment points on degenerate positions:
  newOrder=[]
  for i in range(len(mols)):
    nDegen=0
    for idx,chain in sidechains[i]:
      if idx in degenPts:
        nDegen+=1
    newOrder.append((nDegen,i))
  newOrder.sort()
  
  res=[]
  i=0
  # all the mols without chains in degenerate positions:
  while i<len(newOrder) and newOrder[i][0]==0:
    molIdx = newOrder[i][1]
    tmp = list(sidechains[molIdx])
    tmp.sort()
    res.append((molIdx,tuple(tmp)))
    i+=1

  # all the mols with single chains in degenerate positions:
  chainLocCount=defaultdict(lambda:defaultdict(int))
  posChainLists=defaultdict(list)
  while i<len(newOrder) and newOrder[i][0]==1:
    molIdx = newOrder[i][1]
    # assign the degenerate positions the lowest index:
    tChains=[]
    idx=None
    for tIdx,tChain in sidechains[molIdx]:
      if tIdx not in degenPts:
        tChains.append((tIdx,tChain))
      else:
        idx,chain = tIdx,tChain
    assert idx is not None
    minId=idx
    for matchId,match in enumerate(matches):
      if match[idx]<minId:
        minId=match[idx]
    tChains.append((minId,chain))
    tChains.sort()
    res.append((molIdx,tuple(tChains)))
    chainLocCount[chain][minId]+=1
    posChainLists[minId].append(chain)
    i+=1

  # the remaining molecules have more than one chain in a degenerate position
  while i<len(newOrder) :
    count = newOrder[i][0]
    while i<len(newOrder) and newOrder[i][0]==count:
      molIdx = newOrder[i][1]

      localStates=[]
      seen = set()
      for matchId,match in enumerate(matches):
        chainState = []
        for idx,chain in sidechains[molIdx]:
          pIdx = match[idx]
          chainState.append((idx in degenPts,-chainLocCount[chain][pIdx],pIdx,chain))
        chainState.sort()
        localStates.append(chainState)
      localStates.sort()

      winner = localStates[0]
      tChains=[]
      for a,b,idx,chain in winner:
        tChains.append((idx,chain))
        chainLocCount[chain][idx]+=1
        posChainLists[idx].append(chain)
      tChains.sort()
      res.append((molIdx,tuple(tChains)))
      i+=1
    
  res.sort()
  res = [y for x,y in res]
  
  return res
*/


static
void GetSidechains(const std::vector<ROMOL_SPTR> &mols, const std::vector<ROMOL_SPTR> &cores, 
                   const RGoupDecompositionOptions &options,  
                   std::vector< std::vector< sidechain_t > > &res) {
    
    for(size_t molIdx=0; molIdx < mols.size(); molIdx++) {
        if(NULL == mols[molIdx]) { //never
            /*
//TODO:                          
                  localRes =[(None,'bad molecule')]*len(cores)
                  res.append(tuple(localRes))
            */
            continue;
        }
        const ROMol& mol = *mols[molIdx];
        bool molMatched = false;
    localRes =[]
        for(size_t i=0; i < cores.size(); i++) {
            if(NULL == cores[i]) {
//TODO:                          
                //localRes.append((None,'bad core'));
                continue;
            }
            const ROMol& core = *cores[i];

            std::vector< .... > coreRes;
            std::vector<MatchVectType> tmatches;
            SubstructMatch(mol, core, tmatches);
            if(tmatches.empty()) {
                if(options.Verbose)
                    std::cout<<"molecule "<<molIdx" did not match core "<<i<<"\n";
//TODO:                          
//                coreRes=[(None,'no core matches')]
            }
            else if(tmatches.size() > 1) {
                if(options.Verbose)
                    std::cout<<"core "<<i<<" matches molecule multiple times\n";
                coreRes=[(None,'core matches multiple times')]
            }
            else {
                std::auto_ptr<ROMol> tMol = replaceCore(*mol, *core, true, true);
                if(tMol->getNumAtoms() > 0){
                    molMatched = true;
                    //Fragment tMol:    //tSmi = Chem.MolToSmiles(tMol,True)   tSmi = tSmi.split('.') :
                    std::vector< std::vector<int> > frags;
                    unsigned nFrags = MolOps::getMolFrags(*tMol, frags);
                    for(size_t fi=0; fi < frags.size(); fi++) { // for elem in tSmi:
                        std::map<unsigned, unsigned> newAtomMap;
                        unsigned matches = 0;    // number of attachment points
                        unsigned idx = 0;
                        ROMol* sideChain = new ROMol(); // == cSmi / elem
                        for(size_t ai=0; ai < frags[i].size(); ai++) {
                            const Atom* a = tMol->getAtomWithIdx(frags[fi][ai]);
                            if(0 == a->getElementNum() && LABEL.......... ) // attachment point like [1*] smarts
                            {
                                matches++;  // matches = dummyExpr.findall(elem)   dummyExpr=re.compile(r'\[([0-9]*)\*\]')
                                if( 1 == matches) // first attachment point
                                    idx = frags[fi][ai]; // or ai ???????????????    // ==  idx = matches[0]
                            }
                            else {  // add atom to sidechain (all atoms except attachment points)
                                    // elem = dummyExpr.sub('*',elem)    // REMOVE Labled attachment point [*] ???????
// Is it right ????????????
                                newAtomMap[frags[fi][ai]] = sideChain->addAtom(a->copy(), true, true);
                            }
                        }
                        //add bonds between added atoms
                        std::map<unsigned, unsigned> visitedBonds;
                        for(size_t ai=0; ai < frags[i].size(); ai++) {
                            Atom* a = tMol->getAtomWithIdx(frags[fi][ai]);
                            ROMol::OEDGE_ITER beg,end;
                            for(boost::tie(beg,end) = tMol->getAtomBonds(a); beg!=end; ++beg){
                                const BOND_SPTR bond = (*tMol)[*beg];
                                if(newAtomMap.end() == newAtomMap.find(bond->getBeginAtomIdx())
                                || newAtomMap.end() == newAtomMap.find(bond->getEndAtomIdx())
                                || visitedBonds.end() != visitedBonds.find(bond->getIdx()) )
                                    continue;
                                unsigned ai1 = newAtomMap.at(bond->getBeginAtomIdx());
                                unsigned ai2 = newAtomMap.at(bond->getEndAtomIdx());
                                unsigned bi  = sideChain->addBond(ai1, ai2, bond->getBondType());
                                visitedBonds[bond->getIdx()] = bi;
                            }
                        }
                        
                        if(0 == matches)
                            if(options.Verbose)
                                std::cout<<"sidechain with no attachment point label: "<<elem<<"\n";
                        else { 
                            if(matches > 1)
                                if(options.RejectDoubleAttachments) {
//TODO:                          
                                    //coreRes=((None,'multiple attachment points in sidechain'),)
                                    break;
                                }
                                else
                                    if(options.Verbose)
                                        std::cout<<"sidechain with multiple attachment point labels: "<<elem<<"\n";
/* SKIP:  ???????????????????                            
                          try:
                            /// convert  sideChain to sideChain by AvalonTools  --- SKIP for binary ROMol
                            cSmi = pyAvalonTools.GetCanonSmiles(elem,True)
                          except:
                            pass
                          else:
                            elem = cSmi
*/
//TODO:                          
//                        coreRes.append((idx, ROMOL_SPTR(sideChain)); // attachment point index and fragment/core ?????????????
                        }
                    }                  
                }
            }  
//TODO:                          
            //localRes.append(tuple(coreRes))
        }    
        if( ! molMatched && options.Verbose)
            std::cout<<"molecule "<<molIdx<<" did not match any cores\n";
    }
//TODO:                          
    //res.append(tuple(localRes))
}

static
void ProcessCoreLabels(const std::vector<ROMOL_SPTR> &cores, const RGoupDecompositionOptions &options, 
                             std::vector<ROMOL_SPTR> &resultCores) {

    for(size_t coreIdx=0; coreIdx < cores.size(); coreIdx++) {
        std::vector<unsigned> rvm;
        const ROMol& core = *cores[coreIdx];
        const size_t coreNumAtoms = core.getNumAtoms();
        for(size_t i=0; i < coreNumAtoms; i++)
        {
            const Atom* coreAtom = core.getAtomWithIdx(i);
            if(0 == coreAtom->getAtomicNum()) { // == core.GetSubstructMatches(Chem.MolFromSmiles('[*]'))
                //# use the isotope we read in to set the R label.
                //# this is appropriate because the mol file parser sets R group
                //# labels the same way.
                //# if we come in as SMARTS we use the atom mapping facility
                int label=-1;
                if (coreAtom->HasProp(common_properties::molAtomMapNumber):
                    coreAtom->getProp(common_properties::molAtomMapNumber, label);
                else
                    label = coreAtom->GetIsotope();

                if(0 == coreAtom->getDegree())  // attachment point has no neighbors
                    continue;
        //      if len(nbrs)>1: //        logger.warning('attachment point has more than one neighbor')
                ROMol::ADJ_ITER begin,end;
                boost::tie(begin,end) = core->getAtomNeighbors(atom);  //nbrs = dummy.GetNeighbors()
                for( ; begin!=end; begin++) { //nbr = nbrs[0]
                    //coreAtomNeighbors.push_back(core->getAtomWithIdx(*begin));
                    char str[64];
                    snprintf(str, sizeof(str), "%d", label);
                    core->getAtomWithIdx(*begin)->SetProp('_RLabel',str); //nbr
        //      labelDict[nbr.GetIdx()] = label
                    rmv.push_back(*begin);
                    break; // ??? nbr = nbrs[0]
                }
            }
                
        }

        RWMol* em = new RWMol(core);
        std::sort(rmv.begin(), rmv.end(), std::greater<unsigned>());  // reverse sort
        for(std::vector<unsigned>::const_iterator ri = rmv.begin(); ri != rmv.end(); ri++)
          em.RemoveAtom(*ri);
        resultCores.push_back(ROMOL_SPTR(em));
        //em._labelDict=labelDict
    }
}



//=====================================================================
// Public API implementation:
//=====================================================================

void RGroupDecomposite(const std::vector<ROMOL_SPTR> &mols, const std::vector<ROMOL_SPTR> &srcCores, const RGoupDecompositionOptions &options, 
                             std::vector<ROMOL_SPTR> &results)
{

    std::vector<ROMOL_SPTR>  processedCores;
    std::vector<ROMOL_SPTR>* cores = NULL;

    std::vector<MatchVectType> coreCoreMatches(cores.size());  /// ??? 
    if(options.LabelledCores) {
    ProcessCoreLabels(cores, options, processedCores);
    cores = &processedCores;
    std::vector<MatchVectType> matches;
    for(size_t i=0; i < cores.size(); i++) {
        GetSubstructMatch(*cores[i], *cores[i], matches, false);
        if(matches.empty()) // # pathology happens at times with query features in the core
            for(int ai=0; ai < cores[i].getNumAtoms(); ai++)
                coreCoreMatches[i].push_back(std::pair<int,int>(ai,0)); // ??? matches=(tuple(range(core.GetNumAtoms())),)
        else
            for(unsigned j=0; j < matches.size(); j++)
                for(unsigned k=0; k < matches[j].size(); k++)
                    coreCoreMatches[i].push_back(matches[j][k]); // append
        }
    }
    else
        cores = &srcCores;

    std::vector< std::vector< sidechain_t > > sidechains;  // res: set of Sidechains for each core
    GetSidechains(mols, cores, options, sidechains);
    
    if(options.Symmetrize) {
        if(options.Verbose)
            std::cout<<"Symmetrizing R groups\n";
  
        for(size_t coreIdx=0; coreIdx < cores.size(); coreIdx++) {
            std::vector< sidechain_t > sChains;
            SymmetrizeSidechains(mols, cores[coreIdx], sidechains[coreIdx], options, sChains);

    //// ????????????????????????????????????????????????????????????????? :
            for(size_t i=0; i < sidechains.size(); i++) {  //for i,row in enumerate(res):
                row = list(row)
                row[coreIdx] = sChains[i]
                coreSidechains[i] = row
            }
        }
    }
  
    if(options.LabelledCores) {
        if(options.Verbose)
            std::cout<<"Renumbering to match input ordering\n";
        for(size_t coreIdx=0; coreIdx < cores.size(); coreIdx++) {
            const ROMol& core = *cores[coreIdx];
//      # for symmetric cores where we're labelling, we have to deal with the real PITA
//      # of a core like *C1CC1. Since the dummy gets subsumed into its attached atom
//      # when we do the replacement, this could match C1CC1F three different ways. The one
//      # we get isn't determined by logic, but by atom ordering. We want to present
//      # a sensible result, so we have to do some work:

        const MatchVectType& ccM = coreCoreMatches[coreIdx];
        sidechain_t chains(res[coreIdx].size())
        for(size_t i=0; i < chains.size(); i++)   // chains = [x[coreIdx] for x in res]
            chains[i] = res[coreIdx][i];
/*
      for i,row in enumerate(res):
        possibleMappings=[0]*len(ccM)
        chain = chains[i]
        chain = list(chain)
        for cidx,(idx,smi) in enumerate(chain):
          if idx is None:
            res[i] = (idx,smi)
            chain=None
            break
          for midx,mapping in enumerate(ccM):
            if core.GetAtomWithIdx(mapping[idx]).HasProp('_RLabel'):
              #print '    ',cidx,midx,idx,'->',mapping[idx]
              possibleMappings[midx] += 1
        if chain is not None:
//          #print 'POSS:',possibleMappings
          best = sorted(zip(possibleMappings,range(len(possibleMappings))))
          #raise ValueError(best)
          best = best[-1][1]
          bestMapping = ccM[best]
//          #print 'BEST:',bestMapping

          for cidx,(idx,smi) in enumerate(chain):
            if not core.GetAtomWithIdx(bestMapping[idx]).HasProp('_RLabel'):
              if options.requireLabels or options.nonLabelledSubstituentHandling=='FAIL':
                if not options or not options.silent:
                  logger.warning('skipping molecule %d because it has a substitution at an unmarked position'%i)
                res[i] = (None,'substitution at unmarked position')
                chain=None
                break
              elif options.nonLabelledSubstituentHandling=='IGNORE':
                idx=-1
              else:
                idx+=100
            else:
              idx=int(core.GetAtomWithIdx(bestMapping[idx]).GetProp('_RLabel'))
            if idx>=0:
              chain[cidx]=(idx,smi)
            else:
              chain[cidx]=None
        row = list(row)
        if chain is not None:
          chain = tuple([x for x in chain if x is not None])
          row[coreIdx]=chain
          res[i] = row
*/
        }
//// return mols,res
    }  
}

 
} // end of namespace RDKit




