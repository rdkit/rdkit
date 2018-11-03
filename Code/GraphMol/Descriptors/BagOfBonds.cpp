//
//  Copyright (c) 2018, Guillaume GODIN
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "BagOfBonds.h"
#include "MolData3Ddescriptors.h"
#include <Eigen/Dense>

using namespace Eigen;


struct BagMatrix {
  VectorXd UCM;
  std::vector<std::string > BagTag;
};

namespace RDKit {
namespace Descriptors {
namespace {


VectorXd getUpperMat(MatrixXd mat) {
  VectorXd res(mat.rows()*(mat.cols()+1)/2);
  Index size = mat.rows();
  Index offset = 0;
  for(Index j=0; j<mat.cols(); ++j) {
      res.segment(offset,size) = mat.col(j).tail(size);
      offset += size;
      size--;
  }
  return res;
}

std::list<std::pair<std::string, unsigned int> > EnumeratesAtomsPair(const RDKit::ROMol &mol, 
    bool separateIsotopes, bool abbreviateHIsotopes) {
  std::map<std::string, unsigned int> counts;
  unsigned int nHs = 0;
  std::string key;
  const RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
  for (RDKit::ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    int atNum = (*atomIt)->getAtomicNum();
    key = table->getElementSymbol(atNum);
    if (separateIsotopes) {
      unsigned int isotope = (*atomIt)->getIsotope();
      if (abbreviateHIsotopes && atNum == 1 && (isotope == 2 || isotope == 3)) {
        if (isotope == 2)
          key = "D";
        else
          key = "T";
      }
    }
    if (counts.find(key) != counts.end()) {
      counts[key] += 1;
    } else {
      counts[key] = 1;
    }
    nHs += (*atomIt)->getTotalNumHs();
  }

  if (nHs) {
    key = "H";
    if (counts.find(key) != counts.end()) {
      counts[key] += nHs;
    } else {
      counts[key] = nHs;
    }
  }

  // create the Diagonal single atoms count
  std::list<std::pair<std::string, unsigned int> > ks;
  for (std::map<std::string, unsigned int>::const_iterator countIt = counts.begin();
       countIt != counts.end(); ++countIt) {
    // store the diagonal (individual atoms)
    std::pair<std::string, unsigned int > key = std::make_pair(countIt->first, countIt->second);
    ks.push_back(key);
  }

  unsigned int nbElements = counts.size();
  // enumerate Heterogenous Atom Pairs (simplely number product)
  // heterogens pairs
  for (unsigned int i = 0; i < nbElements - 1 ; i++) {
      std::map<std::string, unsigned int>::const_iterator it1 = counts.begin();
      std::advance(it1,i);
      std::pair<std::string, unsigned int> Values1 = *it1;
      for (unsigned int j = i+1; j < nbElements ; j++) {
        std::map<std::string, unsigned int>::const_iterator it2= counts.begin();
        std::advance(it2,j);
        std::pair<std::string, unsigned int> Values2 = *it2;
        std::string PairString;
        // sorting string order
        if (Values1.first > Values2.first) {
           PairString = Values2.first+Values1.first;
        }
        else {
            PairString = Values1.first+Values2.first;
        }
        unsigned int nbElementPair = Values1.second*Values2.second;
        std::pair<std::string,unsigned int > key = std::make_pair(PairString, nbElementPair);
        ks.push_back(key);
      }
  }

  for (std::map<std::string, unsigned int>::const_iterator countIt = counts.begin();
       countIt != counts.end(); ++countIt) {
    //"homologues":
    unsigned int nbElementAtom = countIt->second;
    unsigned int res =nbElementAtom*(nbElementAtom-1)/2;
    std::string PairString = countIt->first+countIt->first;
    if (res>0) {
      std::pair<std::string,unsigned int > key = std::make_pair(PairString, res);
      ks.push_back(key);
    }
  }
  return ks;

}


 BagMatrix getCoulombMatBags(Eigen::VectorXd numbers, int numatoms, Eigen::MatrixXd Distance3D,  std::vector<std::string > S, int alpha) {
  // 3D distance matrix
  Eigen::MatrixXd ProdV = numbers*numbers.transpose(); // outer products of vector ie (z[i] * z[j])

  if (alpha != 1) {
    Distance3D = Distance3D.array().pow(alpha); // power distance alpha
  }

  Eigen::MatrixXd MyVal = ProdV.cwiseQuotient(Distance3D); // ratio top / dist
  MyVal.diagonal() = 0.5 * numbers.array().pow(2.4); // set the diagonal fix values 
  // get the BagsCode
  std::vector<std::string > BagsCodes;
  std::list<std::pair<std::string, double > > Res;
  std::string s1, s2, key;
  for (unsigned int i = 0; i< numatoms; i++) {
    s1 = S[i];
    BagsCodes.push_back(s1); // adding the diagonal once!

    for (unsigned int j = i+1; j< numatoms; j++ ){
      s2 = S[j];
      if (s1 > s2) {
        key = s2+s1;
      }
      else {
        key = s1+s2;
      }
      BagsCodes.push_back(key);
    }
  }    

  Eigen::VectorXd RES = getUpperMat(MyVal);
  return BagMatrix{RES , BagsCodes}; 
}

///// code to be expose to python directly not simple
// we can parallelized cause we need the max value per key
 std::map<std::string, unsigned int> getBoBSizefromSmiles(std::vector<std::string> smiles) {
    //unsigned int idx = 0;
    //std::vector< std::list < std::pair<std::string, unsigned int> > > myData(smiles.size());
    std::list<std::pair<std::string, unsigned int> > fulldata;
    std::map<std::string, unsigned int> Global;
    std::string key;

    //while (idx < smiles.size()) {
    for (auto const& smi :  smiles) {
      RDKit::ROMol *mol;
      mol = RDKit::SmilesToMol(smi);
      std::list<std::pair<std::string, unsigned int>> data = EnumeratesAtomsPair(*mol, false, false);

      for (auto const& pair :  data) {
          key = pair.first;
          // don't need to find the value just if exists! count instead of find is faster
          if (Global.count(key) > 0) {
            if (pair.second > Global[key]) {
              Global[key] = pair.second;
            }   
          }         
          else {
              Global[key] = pair.second;
          }
      }

      //++idx;
      delete mol;
    }


    return Global;
}

void getBagOfBonds(const RDKit::ROMol &mol, std::vector<double> &res, double *dist3D, 
unsigned int numAtoms, std::map<std::string, unsigned int> MaxBags, int alpha) {

  // initialize the variables for the getCoulombMat
  int numatoms= mol.getNumAtoms();
  double *z = new double[numatoms];
  std::vector<std::string > S;
  for (int i=0; i< numatoms; i++){
      z[i] = mol.getAtomWithIdx(i)->getAtomicNum();
      S.push_back(mol.getAtomWithIdx(i)->getSymbol());
  }

  Eigen::VectorXd numbers = Map<VectorXd>(z, numatoms); // convert the number array to vector
  
  Eigen::MatrixXd Distance3D = Map<MatrixXd>(dist3D, numatoms, numatoms); // convert the result array to matrix (1 column)

  BagMatrix CMBags = getCoulombMatBags(numbers,  numatoms,  Distance3D, S ,alpha);
  // return the BagTag and UCM values as "upper Diagonal elements" only

  // extract the bag CM values index from BagMatrix into sorted vectors padded using MaxBags order & position index!
  std::string key;
  unsigned int sizemax;
  std::vector<double> val;
  for (std::map<std::string, unsigned int>::const_iterator MBags = MaxBags.begin();
       MBags != MaxBags.end(); ++MBags) {
         // take the Master Code Atom Pair
         key = MBags->first;
         std::vector<double> mybagarray;
         unsigned int i=0;
         // loop other the CMBags Code Atom Pair
         for (const auto& element : CMBags.BagTag){ 
            if (element == key) {
                mybagarray.push_back(CMBags.UCM[i]); 
            }
            i++;
         } 

         sizemax = MBags->second;
         //std::cout << key << ": found :" << mybagarray.size() << ", max : " << sizemax << "\n";

         std::vector<double> b;
         if (sizemax-mybagarray.size() > 0) {
           b.clear();
           b.resize(sizemax-mybagarray.size());
         }
         else {
           b.clear();
           b.resize(0);
         }
         if (mybagarray.size() > 0 )
         {
           // descending order
           std::sort(mybagarray.begin(), mybagarray.end());

           std::reverse(mybagarray.begin(), mybagarray.end());
       
           mybagarray.insert(std::end(mybagarray), std::begin(b), std::end(b));
         }
         else {
           mybagarray = b;
         }
        
        val.insert(std::end(val), std::begin(mybagarray), std::end(mybagarray));
  }
    res = val;
}


void getBoBVector(const ROMol &mol, std::vector<double> &res, int confId, unsigned int numAtoms, int alpha,
               std::map<std::string, unsigned int> MaxBags) {
    // 3D distance matrix
    double *dist3D = MolOps::get3DDistanceMat(mol, confId, false, true);
    res.clear();
    res.resize(1);
    getBagOfBonds(mol, res, dist3D, numAtoms,  MaxBags, alpha);
}
}  // end of anonymous namespace

void BagOfBondsVector(const ROMol &mol, std::vector<double> &res, int confId, int alpha,
    std::map<std::string, unsigned int> MaxBags) {
    PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
    unsigned int numAtoms = mol.getNumAtoms();
    // check that here or in the next call ?
    res.clear();
    res.resize(1);
    getBoBVector(mol, res, confId, numAtoms, alpha, MaxBags);
}

std::map<std::string, unsigned int>  BagOfBondsMap(std::vector<std::string>  smiles) {
    std::map<std::string, unsigned int> res =  getBoBSizefromSmiles(smiles);
    return res;
}

}  // namespace Descriptors
}  // namespace RDKit