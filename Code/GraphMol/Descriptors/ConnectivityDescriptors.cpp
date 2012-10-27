// $Id$
//
//  Copyright (C) 2012 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include "MolDescriptors.h"
#include "ConnectivityDescriptors.h"

namespace RDKit{
  namespace Descriptors {
    namespace detail {
      void hkDeltas(const ROMol &mol,std::vector<double> &deltas,bool force){
        PRECONDITION(deltas.size()>=mol.getNumAtoms(),"bad vector size");
        if(!force && mol.hasProp("_connectivityHKDeltas")){
          mol.getProp("_connectivityHKDeltas",deltas);
          return;
        }
        const PeriodicTable *tbl = PeriodicTable::getTable();
        ROMol::VERTEX_ITER atBegin,atEnd;
        boost::tie(atBegin,atEnd) = mol.getVertices();  
        while(atBegin!=atEnd){
          ATOM_SPTR at=mol[*atBegin];
          unsigned int n=at->getAtomicNum();
          if(n<=1){
            deltas[at->getIdx()]=0;
          } else if(n<=10){
            deltas[at->getIdx()]=tbl->getNouterElecs(n)-at->getTotalNumHs();
          } else {
            deltas[at->getIdx()]=double(tbl->getNouterElecs(n)-at->getTotalNumHs())/(n-tbl->getNouterElecs(n)-1);            
          }
          if(deltas[at->getIdx()]!=0.0) deltas[at->getIdx()]=1./sqrt(deltas[at->getIdx()]);
          ++atBegin;
        }
        mol.setProp("_connectivityHKDeltas",deltas);
      }


      void nVals(const ROMol &mol,std::vector<double> &nVs,bool force){
        PRECONDITION(nVs.size()>=mol.getNumAtoms(),"bad vector size");
        if(!force && mol.hasProp("_connectivityNVals")){
          mol.getProp("_connectivityNVals",nVs);
          return;
        }
        const PeriodicTable *tbl = PeriodicTable::getTable();
        ROMol::VERTEX_ITER atBegin,atEnd;
        boost::tie(atBegin,atEnd) = mol.getVertices();  
        while(atBegin!=atEnd){
          ATOM_SPTR at=mol[*atBegin];
          double v=tbl->getNouterElecs(at->getAtomicNum())-at->getTotalNumHs();
          if(v!=0.0){
            v = 1./sqrt(v);
          }
          nVs[at->getIdx()]=v;
          ++atBegin;
        }
        mol.setProp("_connectivityNVals",nVs);
      }


      double chiNv(const ROMol &mol,int n,bool force){
        std::vector<double> hkDs(mol.getNumAtoms());
        detail::hkDeltas(mol,hkDs,force);
        PATH_LIST ps=findAllPathsOfLengthN(mol,n+1,false);
        double res=0.0;
        BOOST_FOREACH(PATH_TYPE p,ps){
          double accum=1.0;
          BOOST_FOREACH(int aidx,p){
            accum*=hkDs[aidx];
          }
          res+=accum;
        }
        return res;
      }
      double chiNn(const ROMol &mol,int n,bool force){
        std::vector<double> nVs(mol.getNumAtoms());
        detail::nVals(mol,nVs,force);
        PATH_LIST ps=findAllPathsOfLengthN(mol,n+1,false);
        double res=0.0;
        BOOST_FOREACH(PATH_TYPE p,ps){
          double accum=1.0;
          BOOST_FOREACH(int aidx,p){
            accum*=nVs[aidx];
          }
          res+=accum;
        }
        return res;
      }

#if 0
      double getAlpha(const Atom &atom,bool &found){
        double res=0.0;
        found=false;
        switch(atom.getAtomicNum()){
        case 1:
          res=1.0;
          found=true;break;
        case 6:
          switch(atom.getHybridization()){
          case Atom::SP:
            res=-0.22;
            found=true;break;
          case Atom::SP2:
            res=-0.13;
            found=true;break;
          case Atom::SP3:
            res=0.00;
            found=true;break;
          default:
            break;
          };
          break;
        case 7:
          switch(atom.getHybridization()){
          case Atom::SP:
            res=-0.29;
            found=true;break;
          case Atom::SP2:
            res=-0.20;
            found=true;break;
          case Atom::SP3:
            res=-0.04;
            found=true;break;
          default:
            break;
          };
          break;
        case 8:
          switch(atom.getHybridization()){
          case Atom::SP2:
            res=-0.20;
            found=true;break;
          case Atom::SP3:
            res=-0.04;
            found=true;break;
          default:
            break;
          };
          break;
        default:
          break;
        }
        case 9:
          switch(atom.getHybridization()){
          case Atom::SP3:
            res=-0.07;
            found=true;break;
          default:
            break;
          };
          break;
        case 15:
          switch(atom.getHybridization()){
          case Atom::SP2:
            res=0.30;
            found=true;break;
          case Atom::SP3:
            res=0.43;
            found=true;break;
          default:
            break;
          };
          break;
        default:
          break;
        }
        case 16:
          switch(atom.getHybridization()){
          case Atom::SP2:
            res=0.22;
            found=true;break;
          case Atom::SP3:
            res=0.35;
            found=true;break;
          default:
            break;
          };
          break;
        default:
          break;
        }
        case 17:
          switch(atom.getHybridization()){
          case Atom::SP3:
            res=0.29;
            found=true;break;
          default:
            break;
          };
          break;
        case 35:
          switch(atom.getHybridization()){
          case Atom::SP3:
            res=0.48;
            found=true;break;
          default:
            break;
          };
          break;
        case 53:
          switch(atom.getHybridization()){
          case Atom::SP3:
            res=0.73;
            found=true;break;
          default:
            break;
          };
          break;
        default:
          break;
        }
        return res;
      }
#endif
    } // end of detail namespace

    double calcChi0v(const ROMol &mol,bool force){
      std::vector<double> hkDs(mol.getNumAtoms());
      detail::hkDeltas(mol,hkDs,force);
      return std::accumulate(hkDs.begin(),hkDs.end(),0.0);
    };
    double calcChi1v(const ROMol &mol,bool force){
      std::vector<double> hkDs(mol.getNumAtoms());
      detail::hkDeltas(mol,hkDs,force);
      
      double res=0.0;
      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
        BOND_SPTR bond = mol[*firstB];
        res += hkDs[bond->getBeginAtomIdx()]*hkDs[bond->getEndAtomIdx()];
        ++firstB;
      }
      return res;
    };
    double calcChi2v(const ROMol &mol,bool force){
      return detail::chiNv(mol,2,force);
    };
    double calcChi3v(const ROMol &mol,bool force){
      return detail::chiNv(mol,3,force);
    };
    double calcChi4v(const ROMol &mol,bool force){
      return detail::chiNv(mol,4,force);
    };

    double calcChi0n(const ROMol &mol,bool force){
      std::vector<double> nVs(mol.getNumAtoms());
      detail::nVals(mol,nVs,force);
      return std::accumulate(nVs.begin(),nVs.end(),0.0);
    };
    double calcChi1n(const ROMol &mol,bool force){
      std::vector<double> nVs(mol.getNumAtoms());
      detail::hkDeltas(mol,nVs,force);
      
      double res=0.0;
      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
        BOND_SPTR bond = mol[*firstB];
        res += nVs[bond->getBeginAtomIdx()]*nVs[bond->getEndAtomIdx()];
        ++firstB;
      }
      return res;
    };
    double calcChi2n(const ROMol &mol,bool force){
      return detail::chiNn(mol,2,force);
    };
    double calcChi3n(const ROMol &mol,bool force){
      return detail::chiNn(mol,3,force);
    };
    double calcChi4n(const ROMol &mol,bool force){
      return detail::chiNn(mol,4,force);
    };

#if 0
    double calcHallKierAlpha(const ROMol &mol){
      const PeriodicTable *tbl = PeriodicTable::getTable();
      double alphaSum=0.0;
      double rC=tbl->getRb0(6);

        ROMol::VERTEX_ITER atBegin,atEnd;
        boost::tie(atBegin,atEnd) = mol.getVertices();  
        while(atBegin!=atEnd){
          ATOM_SPTR at=mol[*atBegin];
          ++atBegin;
          unsigned int n = at->getAtomicNum();
          if(!n) continue;
          bool found;
          double alpha=detail::getAlpha(*(at.get()),found);
          if(!found){
            double rA=tbl->getRb0(n);
            alpha = rA/rC-1.0;
          }
          
        }


      /*
  for atom in m.GetAtoms():
    atNum=atom.GetAtomicNum()
    if not atNum: continue
    symb = atom.GetSymbol()
    alphaV = PeriodicTable.hallKierAlphas.get(symb,None)
    if alphaV is not None:
      hyb = atom.GetHybridization()-2
      if(hyb<len(alphaV)):
        alpha = alphaV[hyb]
        if alpha is None:
          alpha = alphaV[-1]
      else:
        alpha = alphaV[-1]
    else:
      rA = PeriodicTable.nameTable[symb][5]
      alpha = rA/rC - 1
    alphaSum += alpha  
  return alphaSum    
      */      
    };

#endif
    
  } // end of namespace Descriptors
}
