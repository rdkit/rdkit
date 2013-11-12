// $Id$
//
//  Copyright (C) 2004-2012 Greg Landrum and  Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "MolDescriptors.h"
#include "Crippen.h"
#include <iostream>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <boost/lexical_cast.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace RDKit{
  namespace Descriptors {
    extern const std::string defaultParamData;

    void getCrippenAtomContribs(const ROMol &mol,
				std::vector< double > &logpContribs,
				std::vector< double > &mrContribs,
				bool force,
                                std::vector<unsigned int> *atomTypes,
                                std::vector<std::string> *atomTypeLabels
                                ){
      PRECONDITION(logpContribs.size()==mol.getNumAtoms() &&
		   mrContribs.size()==mol.getNumAtoms(),
		   "bad result vector size");
      PRECONDITION((!atomTypes || atomTypes->size()==mol.getNumAtoms()),
                   "bad atomTypes vector");
      PRECONDITION((!atomTypeLabels || atomTypeLabels->size()==mol.getNumAtoms()),
                   "bad atomTypeLabels vector");
      if(!force && mol.hasProp("_crippenLogPContribs")){
	std::vector<double> tmpVect1,tmpVect2;
	mol.getProp("_crippenLogPContribs",tmpVect1);
	mol.getProp("_crippenMRContribs",tmpVect2);
	if(tmpVect1.size()==mol.getNumAtoms() &&
	   tmpVect2.size()==mol.getNumAtoms() ){
	  logpContribs=tmpVect1;
	  mrContribs=tmpVect2;
	  return;
	}
      }
      
      boost::dynamic_bitset<> atomNeeded(mol.getNumAtoms());
      atomNeeded.set();
      const CrippenParamCollection *params=CrippenParamCollection::getParams();
      for(CrippenParamCollection::ParamsVect::const_iterator it=params->begin();
	  it!=params->end(); ++it){
	std::vector<MatchVectType> matches;
	SubstructMatch(mol,*(it->dp_pattern.get()),matches,
		       false,true);
	for(std::vector<MatchVectType>::const_iterator matchIt=matches.begin();
	    matchIt!=matches.end();++matchIt){
	  int idx=(*matchIt)[0].second;
	  if(atomNeeded[idx]){
	    atomNeeded[idx]=0;
	    logpContribs[idx] = it->logp;
	    mrContribs[idx] = it->mr;
            if(atomTypes) (*atomTypes)[idx]=it->idx;
            if(atomTypeLabels) (*atomTypeLabels)[idx]=it->label;
	  }
	}
	// no need to keep matching stuff if we already found all the atoms:
	if(atomNeeded.none()) break;
      }
      mol.setProp("_crippenLogPContribs",logpContribs,true);
      mol.setProp("_crippenMRContribs",mrContribs,true);
    }
    void calcCrippenDescriptors(const ROMol &mol,double &logp,double &mr,bool includeHs,
				bool force){
      if(!force && mol.hasProp("_crippenLogP")){
	mol.getProp("_crippenLogP",logp);
	mol.getProp("_crippenMR",mr);
	return;
      }

      // this isn't as bad as it looks, we aren't actually going
      // to harm the molecule in any way!
      ROMol *workMol=const_cast<ROMol *>(&mol);
      if(includeHs){
	workMol = MolOps::addHs(mol,false,false);
      }
      std::vector<double> logpContribs(workMol->getNumAtoms());
      std::vector<double> mrContribs(workMol->getNumAtoms());
      getCrippenAtomContribs(*workMol,logpContribs,mrContribs,force);
      logp=0.0;
      for(std::vector<double>::const_iterator iter=logpContribs.begin();
	  iter!=logpContribs.end();++iter){
	logp+= *iter;
      }
      mr=0.0;
      for(std::vector<double>::const_iterator iter=mrContribs.begin();
	  iter!=mrContribs.end();++iter){
	mr+= *iter;
      }

      if(includeHs){
	delete workMol;
      }

      mol.setProp("_crippenLogP",logp,true);
      mol.setProp("_crippenMR",mr,true);
    };

    typedef boost::flyweight<boost::flyweights::key_value<std::string,CrippenParamCollection>,
                             boost::flyweights::no_tracking > param_flyweight;

    const CrippenParamCollection *CrippenParamCollection::getParams(const std::string &paramData){
      const CrippenParamCollection *res = &(param_flyweight(paramData).get());
      return res;
    }
    CrippenParamCollection::CrippenParamCollection(const std::string &paramData){
      std::string params;
      boost::char_separator<char> tabSep("\t","",boost::keep_empty_tokens);
      if(paramData=="") params=defaultParamData;
      else params=paramData;
      std::istringstream inStream(params);

      std::string inLine=RDKit::getLine(inStream);
      unsigned int idx=0;
      while(!inStream.eof()){
	if(inLine[0] != '#'){
	  CrippenParams paramObj;
          paramObj.idx=idx++;
	  tokenizer tokens(inLine,tabSep);
	  tokenizer::iterator token=tokens.begin();
	  
	  paramObj.label=*token;
	  ++token;
	  paramObj.smarts=*token;
	  ++token;
	  if(*token!=""){
	    paramObj.logp=boost::lexical_cast<double>(*token);
	  } else {
	    paramObj.logp=0.0;
	  }
	  ++token;
	  if(*token!=""){
	    try{
	      paramObj.mr=boost::lexical_cast<double>(*token);
	    } catch (boost::bad_lexical_cast){
	      paramObj.mr=0.0;
	    }
	  } else {
	      paramObj.mr=0.0;
	  }
	  paramObj.dp_pattern=boost::shared_ptr<const ROMol>(SmartsToMol(paramObj.smarts));
	  d_params.push_back(paramObj);
	}
	inLine = RDKit::getLine(inStream);
      }
    }
    
    CrippenParams::~CrippenParams() {
      dp_pattern.reset();
    }

const std::string defaultParamData=
"#ID	SMARTS	logP	MR	Notes/Questions\n"
"C1	[CH4]	0.1441	2.503	\n"
"C1	[CH3]C	0.1441	2.503	\n"
"C1	[CH2](C)C	0.1441	2.503	\n"
"C2	[CH](C)(C)C	0	2.433	\n"
"C2	[C](C)(C)(C)C	0	2.433	\n"
"C3	[CH3][N,O,P,S,F,Cl,Br,I]	-0.2035	2.753	\n"
"C3	[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]	-0.2035	2.753	\n"
"C4	[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]	-0.2051	2.731	\n"
"C4	[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]	-0.2051	2.731	\n"
"C5	[C]=[!C;A;!#1]	-0.2783	5.007	\n"
"C6	[CH2]=C	0.1551	3.513	\n"
"C6	[CH1](=C)[A;!#1]	0.1551	3.513	\n"
"C6	[CH0](=C)([A;!#1])[A;!#1]	0.1551	3.513	\n"
"C6	[C](=C)=C	0.1551	3.513	\n"
"C7	[CX2]#[A;!#1]	0.0017	3.888	\n"
"C8	[CH3]c	0.08452	2.464	\n"
"C9	[CH3]a	-0.1444	2.412	\n"
"C10	[CH2X4]a	-0.0516	2.488	\n"
"C11	[CHX4]a	0.1193	2.582	\n"
"C12	[CH0X4]a	-0.0967	2.576	\n"
"C13	[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]	-0.5443	4.041	\n"
"C14	[c][#9]	0	3.257	\n"
"C15	[c][#17]	0.245	3.564	\n"
"C16	[c][#35]	0.198	3.18	\n"
"C17	[c][#53]	0	3.104	\n"
"C18	[cH]	0.1581	3.35	\n"
"C19	[c](:a)(:a):a	0.2955	4.346	\n"
"C20	[c](:a)(:a)-a	0.2713	3.904	\n"
"C21	[c](:a)(:a)-C	0.136	3.509	\n"
"C22	[c](:a)(:a)-N	0.4619	4.067	\n"
"C23	[c](:a)(:a)-O	0.5437	3.853	\n"
"C24	[c](:a)(:a)-S	0.1893	2.673	\n"
"C25	[c](:a)(:a)=[C,N,O]	-0.8186	3.135	\n"
"C26	[C](=C)(a)[A;!#1]	0.264	4.305	\n"
"C26	[C](=C)(c)a	0.264	4.305	\n"
"C26	[CH1](=C)a	0.264	4.305	\n"
"C26	[C]=c	0.264	4.305	\n"
"C27	[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]	0.2148	2.693	\n"
"CS	[#6]	0.08129	3.243	\n"
"H1	[#1][#6,#1]	0.123	1.057	\n"
"H2	[#1]O[CX4,c]	-0.2677	1.395	\n"
"H2	[#1]O[!#6;!#7;!#8;!#16]	-0.2677	1.395	\n"
"H2	[#1][!#6;!#7;!#8]	-0.2677	1.395	\n"
"H3	[#1][#7]	0.2142	0.9627	\n"
"H3	[#1]O[#7]	0.2142	0.9627	\n"
"H4	[#1]OC=[#6,#7,O,S]	0.298	1.805	\n"
"H4	[#1]O[O,S]	0.298	1.805	\n"
"HS	[#1]	0.1125	1.112	\n"
"N1	[NH2+0][A;!#1]	-1.019	2.262	\n"
"N2	[NH+0]([A;!#1])[A;!#1]	-0.7096	2.173	\n"
"N3	[NH2+0]a	-1.027	2.827	\n"
"N4	[NH1+0]([!#1;A,a])a	-0.5188	3	\n"
"N5	[NH+0]=[!#1;A,a]	0.08387	1.757	\n"
"N6	[N+0](=[!#1;A,a])[!#1;A,a]	0.1836	2.428	\n"
"N7	[N+0]([A;!#1])([A;!#1])[A;!#1]	-0.3187	1.839	\n"
"N8	[N+0](a)([!#1;A,a])[A;!#1]	-0.4458	2.819	\n"
"N8	[N+0](a)(a)a	-0.4458	2.819	\n"
"N9	[N+0]#[A;!#1]	0.01508	1.725	\n"
"N10	[NH3,NH2,NH;+,+2,+3]	-1.95		\n"
"N11	[n+0]	-0.3239	2.202	\n"
"N12	[n;+,+2,+3]	-1.119		\n"
"N13	[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]	-0.3396	0.2604	\n"
"N13	[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]	-0.3396	0.2604	\n"
"N13	[NH0;+,+2,+3](=[#6])=[#7]	-0.3396	0.2604	\n"
"N14	[N;+,+2,+3]#[A;!#1]	0.2887	3.359	\n"
"N14	[N;-,-2,-3]	0.2887	3.359	\n"
"N14	[N;+,+2,+3](=[N;-,-2,-3])=N	0.2887	3.359	\n"
"NS	[#7]	-0.4806	2.134	\n"
"O1	[o]	0.1552	1.08	\n"
"O2	[OH,OH2]	-0.2893	0.8238	\n"
"O3	[O]([A;!#1])[A;!#1]	-0.0684	1.085	\n"
"O4	[O](a)[!#1;A,a]	-0.4195	1.182	\n"
"O5	[O]=[#7,#8]	0.0335	3.367	\n"
"O5	[OX1;-,-2,-3][#7]	0.0335	3.367	\n"
"O6	[OX1;-,-2,-2][#16]	-0.3339	0.7774	\n"
"O6	[O;-0]=[#16;-0]	-0.3339	0.7774	\n"
"O12	[O-]C(=O)	-1.326		\"order flip here intentional\"\n"
"O7	[OX1;-,-2,-3][!#1;!N;!S]	-1.189	0	\n"
"O8	[O]=c	0.1788	3.135	\n"
"O9	[O]=[CH]C	-0.1526	0	\n"
"O9	[O]=C(C)([A;!#1])	-0.1526	0	\n"
"O9	[O]=[CH][N,O]	-0.1526	0	\n"
"O9	[O]=[CH2]	-0.1526	0	\n"
"O9	[O]=[CX2]=O	-0.1526	0	\n"
"O10	[O]=[CH]c	0.1129	0.2215	\n"
"O10	[O]=C([C,c])[a;!#1]	0.1129	0.2215	\n"
"O10	[O]=C(c)[A;!#1]	0.1129	0.2215	\n"
"O11	[O]=C([!#1;!#6])[!#1;!#6]	0.4833	0.389	\n"
"OS	[#8]	-0.1188	0.6865	\n"
"F	[#9-0]	0.4202	1.108	\n"
"Cl	[#17-0]	0.6895	5.853	\n"
"Br	[#35-0]	0.8456	8.927	\n"
"I	[#53-0]	0.8857	14.02	\n"
"Hal	[#9,#17,#35,#53;-]	-2.996		\n"
"Hal	[#53;+,+2,+3]	-2.996		\n"
"Hal	[+;#3,#11,#19,#37,#55]	-2.996		\"Footnote h indicates these should be here?\"\n"
"P	[#15]	0.8612	6.92	\n"
"S2	[S;-,-2,-3,-4,+1,+2,+3,+5,+6]	-0.0024	7.365	\"Order flip here is intentional\"\n"
"S2	[S-0]=[N,O,P,S]	-0.0024	7.365	\"Expanded definition of (pseudo-)ionic S\"\n"
"S1	[S;A]	0.6482	7.591	\"Order flip here is intentional\"\n"
"S3	[s;a]	0.6237	6.691	\n"
"Me1	[#3,#11,#19,#37,#55]	-0.3808	5.754	\n"
"Me1	[#4,#12,#20,#38,#56]	-0.3808	5.754	\n"
"Me1	[#5,#13,#31,#49,#81]	-0.3808	5.754	\n"
"Me1	[#14,#32,#50,#82]	-0.3808	5.754	\n"
"Me1	[#33,#51,#83]	-0.3808	5.754	\n"
"Me1	[#34,#52,#84]	-0.3808	5.754	\n"
"Me2	[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]	-0.0025		\n"
"Me2	[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]	-0.0025		\n"
"Me2	[#72,#73,#74,#75,#76,#77,#78,#79,#80]	-0.0025		\n";

  } // end of namespace Descriptors
}
