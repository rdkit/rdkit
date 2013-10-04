/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%{
#include <GraphMol/Descriptors/MolDescriptors.h>
%}

%include <GraphMol/Descriptors/MolDescriptors.h>
%include <GraphMol/Descriptors/Lipinski.h>
%include <GraphMol/Descriptors/MolSurf.h>
%include <GraphMol/Descriptors/ConnectivityDescriptors.h>
%include <GraphMol/Descriptors/MQN.h>

%inline %{
  std::pair<double,double> calcCrippenDescriptors(const RDKit::ROMol &mol,
                                                   bool includeHs=true,bool force=false) {
    std::pair<double,double> res;
    RDKit::Descriptors::calcCrippenDescriptors(mol, res.first, res.second, includeHs, force);
    return res;
  }
  double calcMolLogP(const RDKit::ROMol &mol){
    double logp,mr;
    RDKit::Descriptors::calcCrippenDescriptors(mol,logp,mr);
    return logp;
  }
  double calcMolMR(const RDKit::ROMol &mol){
    double logp,mr;
    RDKit::Descriptors::calcCrippenDescriptors(mol,logp,mr);
    return mr;
  }

  std::vector<double> calcUSR(const RDKit::ROMol &mol,int confId=-1){
    std::vector<double> res(12,0);
    RDKit::Descriptors::USR(mol,res,confId);
    return res;
  }
  std::vector<double> calcUSRCAT(const RDKit::ROMol &mol,int confId=-1){
    std::vector<double> res(60,0);
    std::vector<std::vector<unsigned int> > atomIds(0);
    RDKit::Descriptors::USRCAT(mol,res,atomIds,confId);
    return res;
  }
  double calcUSRScore(const std::vector<double> &v1,const std::vector<double> &v2,
                      const std::vector<double> &weights){
    return RDKit::Descriptors::calcUSRScore(v1,v2,weights);
  }
  double calcUSRScore(const std::vector<double> &v1,const std::vector<double> &v2){
    const std::vector<double> weights(v1.size()/12,1.0);
    return RDKit::Descriptors::calcUSRScore(v1,v2,weights);
  }
%}


