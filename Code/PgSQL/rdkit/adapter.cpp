// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
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
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>
#include <boost/integer_traits.hpp>

#include "rdkit.h"

using namespace std;
using namespace RDKit;

const unsigned int SSS_FP_SIZE=1024;
const unsigned int LAYERED_FP_SIZE=1024;
const unsigned int MORGAN_FP_SIZE=1024;
const unsigned int HASHED_PAIR_FP_SIZE=2048;
class ByteA : public std::string {
	public:
		ByteA() : string() {};
		ByteA(bytea *b) : string(VARDATA(b), VARSIZE(b)-VARHDRSZ) {};
		ByteA(string& s) : string(s) {};

		/*
 		 * Convert string to bytea. Convertaion is in pgsql's memory
		 */	
		bytea*  toByteA() {
			bytea *res;
			int len;

			len = this->size();
     			res = (bytea*)palloc( VARHDRSZ + len );
        		memcpy(VARDATA(res), this->data(), len);
        		SET_VARSIZE(res, VARHDRSZ + len);
			 
			return res;
		};

		/* Just the copy of string's method */
		ByteA& operator=(const string& __str) {return (ByteA&)this->assign(__str);};
};

/*
 * Constant io
 */
static string StringData;

/*
 * Real sparse vector
 */

typedef SparseIntVect<boost::uint32_t> SparseFP;


/*******************************************
 *        ROMol transformation             *
 *******************************************/

extern "C"  void    
freeCROMol(CROMol data) {
	ROMol   *mol = (ROMol*)data;
	delete mol;
}

extern "C" CROMol 
constructROMol(Mol *data) {
	ROMol	*mol = new ROMol();
	
	try {
		ByteA b(data);
		MolPickler::molFromPickle(b, mol);
	} catch (MolPicklerException& e) {
		elog(ERROR, "molFromPickle: %s", e.message());
	} catch (...) {
		elog(ERROR, "constructROMol: Unknown exception");
	}
	
	return (CROMol)mol;	
}

extern "C" Mol* 
deconstructROMol(CROMol data) {
	ROMol   *mol = (ROMol*)data;
	ByteA	b;

	try {
		MolPickler::pickleMol(mol, b);
	} catch (MolPicklerException& e) {
		elog(ERROR, "pickleMol: %s", e.message());
	} catch (...) {
		elog(ERROR, "deconstructROMol: Unknown exception");
	}

	return (Mol*)b.toByteA();
}

extern "C" CROMol 
parseMolText(char *data,bool asSmarts) {
  ROMol   *mol = NULL;

  try {
    StringData.assign(data);
    if(!asSmarts){
      mol = SmilesToMol(StringData);
    } else {
      mol = SmartsToMol(StringData);
    }
  } catch (...) {
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("problem generating molecule from smiles '%s'",data)));
  }
  if(mol==NULL){
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("smiles '%s' could not be parsed",data)));
  }

  return (CROMol)mol;
}

extern "C" bool
isValidSmiles(char *data) {
  ROMol   *mol = NULL;
  bool res;
  try {
    StringData.assign(data);
    mol = SmilesToMol(StringData);
  } catch (...) {
    mol=NULL;
  }
  if(mol==NULL){
    res=false;
  } else {
    res=true;
    delete mol;
  }
  return res;
}

extern "C" bool
isValidSmarts(char *data) {
  ROMol   *mol = NULL;
  bool res;
  try {
    StringData.assign(data);
    mol = SmartsToMol(StringData);
  } catch (...) {
    mol=NULL;
  }
  if(mol==NULL){
    res=false;
  } else {
    res=true;
    delete mol;
  }
  return res;
}

extern "C" char *
makeMolText(CROMol data, int *len,bool asSmarts) {
	ROMol   *mol = (ROMol*)data;

	try {
          if(!asSmarts){
            StringData = MolToSmiles(*mol, true);
          } else {
            StringData = MolToSmarts(*mol, false);
          }
	} catch (...) {
		elog(ERROR, "makeMolText: Unknown exception");
	}	

	*len = StringData.size();
	return (char*)StringData.c_str();		
}

extern "C" bytea* 
makeMolSign(CROMol data) {
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;
	bytea			*ret = NULL;

	try {
          res = RDKit::LayeredFingerprintMol(*mol,0x07,1,7,SSS_FP_SIZE);
          ret = makeSignatureBitmapFingerPrint((MolBitmapFingerPrint)res);
          delete res;
          res=0;
	} catch (...) {
          elog(ERROR, "makeMolSign: Unknown exception");
          if(res) delete res;
	}
	
	return ret;
}

extern "C" int
MolSubstruct(CROMol i, CROMol a) {
	ROMol *im = (ROMol*)i;
	ROMol *am = (ROMol*)a;
	RDKit::MatchVectType matchVect;

	return RDKit::SubstructMatch(*im,*am,matchVect); 
}


/*******************************************
 *     Molecule operations                 *
 *******************************************/
extern "C" double
MolAMW(CROMol i){
  return RDKit::Descriptors::CalcAMW(*(ROMol*)i,false);
}
extern "C" double
MolLogP(CROMol i){
  double logp,mr;
  RDKit::Descriptors::CalcCrippenDescriptors(*(ROMol*)i,logp,mr);
  return logp;
}
extern "C" int
MolHBA(CROMol i){
  const ROMol *im = (ROMol*)i;
  int res=0;
  for(ROMol::ConstAtomIterator iter=im->beginAtoms();
      iter!=im->endAtoms();++iter){
    if((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) ++res;
  }
  return res;
}
extern "C" int
MolHBD(CROMol i){
  const ROMol *im = (ROMol*)i;
  int res=0;
  for(ROMol::ConstAtomIterator iter=im->beginAtoms();
      iter!=im->endAtoms();++iter){
    if(((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) && (*iter)->getTotalNumHs()>0) ++res;
  }
  return res;
}
extern "C" int
MolNumAtoms(CROMol i){
  const ROMol *im = (ROMol*)i;
  return im->getNumAtoms(false);
}
extern "C" int
MolNumHeavyAtoms(CROMol i){
  const ROMol *im = (ROMol*)i;
  return im->getNumAtoms(true);
}



/*******************************************
 *     MolBitmapFingerPrint transformation *
 *******************************************/

extern "C"  void    
freeMolBitmapFingerPrint(MolBitmapFingerPrint data) {
	ExplicitBitVect   *fp = (ExplicitBitVect*)data;
	delete fp;
}

extern "C" MolBitmapFingerPrint 
constructMolBitmapFingerPrint(BitmapFingerPrint *data) {
	ExplicitBitVect *ebv=NULL;
	
	try {
 		ebv = new ExplicitBitVect(VARDATA(data), VARSIZE(data) - VARHDRSZ);
	} catch (...) {
		elog(ERROR, "constructMolFingerPrint: Unknown exception");
	}
	
	return (MolBitmapFingerPrint)ebv;	
}

extern "C" BitmapFingerPrint * 
deconstructMolBitmapFingerPrint(MolBitmapFingerPrint data) {
	ExplicitBitVect *ebv = (ExplicitBitVect*)data;
	ByteA		 b;

	try {
		b = ebv->toString();
	} catch (...) {
		elog(ERROR, "deconstructMolFingerPrint: Unknown exception");
	}

	return b.toByteA();
}

extern "C" bytea *
makeSignatureBitmapFingerPrint(MolBitmapFingerPrint data) {
	ExplicitBitVect *ebv = (ExplicitBitVect*)data;
	int	numBits = ebv->getNumBits(),
		i,
		numBytes;
	bytea	*res;
	unsigned char *s;

	numBytes = VARHDRSZ + (numBits/8);
	if ( (numBits % 8) != 0 ) numBytes++;
		
	res = (bytea*)palloc0(numBytes);
	SET_VARSIZE(res, numBytes);
	s = (unsigned char *)VARDATA(res);

	for(i=0; i<numBits; i++)
		if (ebv->getBit(i))
			s[ i/8 ]  |= 1 << (i % 8); 

	return res;	
}

extern "C" int
MolBitmapFingerPrintSize(MolBitmapFingerPrint a) {
	ExplicitBitVect	*f = (ExplicitBitVect*)a;

	return f->getNumBits();
}

extern "C" double
calcBitmapTanimotoSml(MolBitmapFingerPrint a, MolBitmapFingerPrint b) {
	double res=0.0;

	/*
     * Nsame / (Na + Nb - Nsame)
     */
	
	try {
		res = TanimotoSimilarity(*(ExplicitBitVect*)a, *(ExplicitBitVect*)b);
	} catch (ValueErrorException& e) {
		elog(ERROR, "TanimotoSimilarity: %s", e.message().c_str());
	} catch (...) {
		elog(ERROR, "calcBitmapTanimotoSml: Unknown exception");
	}

	return res;
}

extern "C" double
calcBitmapDiceSml(MolBitmapFingerPrint a, MolBitmapFingerPrint b) {
	double res=0.0;

	/*
     * 2 * Nsame / (Na + Nb)
     */
	
	try {
		res = DiceSimilarity(*(ExplicitBitVect*)a, *(ExplicitBitVect*)b);
	} catch (ValueErrorException& e) {
		elog(ERROR, "DiceSimilarity: %s", e.message().c_str());
	} catch (...) {
		elog(ERROR, "calcTanimotoSml: Unknown exception");
	}

	return res;
}


/*******************************************
 *     MolSparseFingerPrint transformation *
 *******************************************/

extern "C"  void    
freeMolSparseFingerPrint(MolSparseFingerPrint data) {
	SparseFP   *fp = (SparseFP*)data;
	delete fp;
}

extern "C" MolSparseFingerPrint 
constructMolSparseFingerPrint(SparseFingerPrint *data) {
	SparseFP *ebv = NULL;
	
	try {
 		ebv = new SparseFP(VARDATA(data), VARSIZE(data) - VARHDRSZ);
	} catch (...) {
		elog(ERROR, "constructMolFingerPrint: Unknown exception");
	}
	
	return (MolSparseFingerPrint)ebv;	
}

extern "C" SparseFingerPrint * 
deconstructMolSparseFingerPrint(MolSparseFingerPrint data) {
	SparseFP *ebv = (SparseFP*)data;
	ByteA		 b;

	try {
		b = ebv->toString();
	} catch (...) {
		elog(ERROR, "deconstructMolFingerPrint: Unknown exception");
	}

	return b.toByteA();
}

extern "C" bytea *
makeSignatureSparseFingerPrint(MolSparseFingerPrint data, int numBits) {
	SparseFP *v = (SparseFP*)data;
	int	n,
		numBytes;
	bytea	*res;
	unsigned char *s;
	SparseFP::StorageType::const_iterator iter;

	numBytes = VARHDRSZ + (numBits/8);
	if ( (numBits % 8) != 0 ) numBytes++;
		
	res = (bytea*)palloc0(numBytes);
	SET_VARSIZE(res, numBytes);
	s = (unsigned char *)VARDATA(res);


	for(iter = v->getNonzeroElements().begin(); iter != v->getNonzeroElements().end(); iter++)
	{
		n = iter->first % numBits;
		s[ n/8 ]  |= 1 << (n % 8);
	}

	return res;	
}

extern "C" bytea * 
makeLowSparseFingerPrint(MolSparseFingerPrint data, int numInts) {
	SparseFP *v = (SparseFP*)data;
	int		numBytes;
	bytea	*res;
	IntRange *s;
	int		n;
	SparseFP::StorageType::const_iterator iter;

	numBytes = VARHDRSZ + (numInts * sizeof(IntRange));
		
	res = (bytea*)palloc0(numBytes);
	SET_VARSIZE(res, numBytes);
	s = (IntRange *)VARDATA(res);


	for(iter = v->getNonzeroElements().begin(); iter != v->getNonzeroElements().end(); iter++)
	{
                uint32 iterV=(uint32)iter->second;
		n = iter->first % numInts;

		if (iterV > INTRANGEMAX){
#if 0
			elog(ERROR, "sparse fingerprint is too big, increase INTRANGEMAX in rdkit.h");
#else
                        iterV=INTRANGEMAX;
#endif
                }
                
		if (s[ n ].low == 0 || s[ n ].low > iterV)
                  s[ n ].low = iterV;
		if (s[ n ].high < iterV)
                  s[ n ].high = iterV;
	}

	return res;	
}

extern "C" void
countOverlapValues(bytea * sign, MolSparseFingerPrint data, int numBits, 
		int * sum, int * overlapSum, int * overlapN) 
{
	SparseFP *v = (SparseFP*)data;
	SparseFP::StorageType::const_iterator iter;

	*sum = *overlapSum = *overlapN = 0;

	if (sign)
	{
		unsigned char *s = (unsigned char *)VARDATA(sign);
		int		n;

		for(iter = v->getNonzeroElements().begin(); iter != v->getNonzeroElements().end(); iter++)
		{
			*sum += iter->second;
			n = iter->first % numBits;
			if ( s[n/8] & (1 << (n % 8)) )
			{
				*overlapSum += iter->second;
				*overlapN += 1;
			}
		}
	}
	else
	{
		/* Assume, sign has only true bits */
		for(iter = v->getNonzeroElements().begin(); iter != v->getNonzeroElements().end(); iter++)
			*sum += iter->second;

		*overlapSum = *sum;
		*overlapN = v->getNonzeroElements().size(); 
	}
}

extern "C" void 
countLowOverlapValues(bytea * sign, MolSparseFingerPrint data, int numInts,
        int * querySum, int *keySum, int * overlapUp, int * overlapDown)
{
	SparseFP *v = (SparseFP*)data;
	SparseFP::StorageType::const_iterator iter;
	IntRange *s = (IntRange *)VARDATA(sign);
	int		n;

	*querySum = *keySum = *overlapUp = *overlapDown = 0;

	for(iter = v->getNonzeroElements().begin(); iter != v->getNonzeroElements().end(); iter++)
	{
		*querySum += iter->second;
		n = iter->first % numInts;
		if (s[n].low == 0) 
		{
			Assert(s[n].high == 0);
			continue;
		}

		*overlapDown += Min(s[n].low, (uint32)iter->second);
		*overlapUp += Min(s[n].high, (uint32)iter->second);
	}

	Assert(*overlapDown <= *overlapUp);

	for(n=0;n<numInts;n++) 
	{
		*keySum += s[n].low;
		if (s[n].low != s[n].high)
			*keySum += s[n].high; /* there is at least two key mapped into current backet */
	}

	Assert(*overlapUp <= *keySum);
}

extern "C" double
calcSparseTanimotoSml(MolSparseFingerPrint a, MolSparseFingerPrint b) {
	double res = -1.0;

	/*
     * Nsame / (Na + Nb - Nsame)
     */
	
	try {
		res = TanimotoSimilarity(*(SparseFP*)a, *(SparseFP*)b);
	} catch (ValueErrorException& e) {
		elog(ERROR, "TanimotoSimilarity: %s", e.message().c_str());
	} catch (...) {
		elog(ERROR, "calcSparseTanimotoSml: Unknown exception");
	}

	return res;
}

extern "C" double
calcSparseDiceSml(MolSparseFingerPrint a, MolSparseFingerPrint b) {
	double res = -1.0;

	/*
     * 2 * Nsame / (Na + Nb)
     */
	
	try {
		res = DiceSimilarity(*(SparseFP*)a, *(SparseFP*)b);
	} catch (ValueErrorException& e) {
		elog(ERROR, "DiceSimilarity: %s", e.message().c_str());
	} catch (...) {
		elog(ERROR, "calcSparseDiceSml: Unknown exception");
	}

	return res;
}

extern "C" double
calcSparseStringDiceSml(const char *a, unsigned int sza, const char *b, unsigned int szb) {
  const unsigned char *t1=(const unsigned char *)a;
  const unsigned char *t2=(const unsigned char *)b;

  boost::uint32_t tmp;
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  if(tmp!=(boost::uint32_t)ci_SPARSEINTVECT_VERSION){
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 1");
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(tmp!=(boost::uint32_t)ci_SPARSEINTVECT_VERSION){
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 2");
  }

  // check the element size:
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  if(tmp!=sizeof(boost::uint32_t)){
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 1 -> uint32_t");
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(tmp!=sizeof(boost::uint32_t)){
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 2 -> uint32_t");
  }
 
  double res=0.;
  // start reading:
  boost::uint32_t len1,len2;
  len1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  len2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(len1!=len2){
    elog(ERROR, "attempt to compare fingerprints of different length");
  }

  boost::uint32_t nElem1,nElem2;
  nElem1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  nElem2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);

  if(!nElem1 || !nElem2){
    return 0.0;
  }

  double v1Sum=0,v2Sum=0,numer=0;
  boost::uint32_t idx1=0;
  boost::int32_t v1;
  boost::uint32_t idx2=0;
  boost::int32_t v2;
  idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
  t1+=sizeof(boost::int32_t);
  nElem1--;
  v1Sum += v1;

  idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
  t2+=sizeof(boost::int32_t);
  nElem2--;
  v2Sum += v2;

  while(1){
    while(nElem2 && idx2<idx1){
      idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
      t2+=sizeof(boost::uint32_t);
      v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
      t2+=sizeof(boost::int32_t);
      nElem2--;
      v2Sum += v2;
    }
    if(idx2==idx1 ){
      //std::cerr<<"   --- "<<idx1<<" "<<v1<<" - "<<idx2<<" "<<v2<<std::endl;
      numer += std::min(v1,v2);
    }
    if(nElem1){
      idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
      t1+=sizeof(boost::uint32_t);
      v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
      t1+=sizeof(boost::int32_t);
      nElem1--;
      v1Sum += v1;
    } else {
      break;
    }
  }
  while(nElem2){
    idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
    t2+=sizeof(boost::uint32_t);
    v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
    t2+=sizeof(boost::int32_t);
    nElem2--;
    v2Sum += v2;
  }
  double denom=v1Sum+v2Sum;
  if(fabs(denom)<1e-6){
    res=0.0;
  } else {
    res = 2.*numer/denom;
  }

  return res;
}


extern "C" MolSparseFingerPrint 
addSFP(MolSparseFingerPrint a, MolSparseFingerPrint b) {
	SparseFP	*res=NULL;
	try {
          SparseFP tmp=(*(SparseFP*)a+*(SparseFP*)b);
          res=(SparseFP*)new SparseFP(tmp);
	} catch (...) {
		elog(ERROR, "addSFP: Unknown exception");
	}
	return (MolSparseFingerPrint)res;
}

extern "C" MolSparseFingerPrint 
subtractSFP(MolSparseFingerPrint a, MolSparseFingerPrint b) {
	SparseFP	*res=NULL;
	try {
          SparseFP tmp=(*(SparseFP*)a-*(SparseFP*)b);
          res=(SparseFP*)new SparseFP(tmp);
	} catch (...) {
		elog(ERROR, "addSFP: Unknown exception");
	}
	return (MolSparseFingerPrint)res;
}



/*
 * Mol -> fp
 */
extern "C" MolBitmapFingerPrint 
makeLayeredBFP(CROMol data) {
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;

	try {
          res = RDKit::LayeredFingerprintMol(*mol,0xFFFFFFFF,1,7,LAYERED_FP_SIZE);
	} catch (...) {
		elog(ERROR, "makeLayeredBFP: Unknown exception");
                if(res) delete res;
                res=NULL;
	}
	
	return (MolBitmapFingerPrint)res;
}

extern "C" MolBitmapFingerPrint 
makeRDKitBFP(CROMol data) {
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;

	try {
          res = RDKit::RDKFingerprintMol(*mol,1,5,LAYERED_FP_SIZE,1);
	} catch (...) {
		elog(ERROR, "makeRDKitBFP: Unknown exception");
                if(res) delete res;
                res=NULL;
	}
	
	return (MolBitmapFingerPrint)res;
}

extern "C" MolSparseFingerPrint 
makeMorganSFP(CROMol data, int radius) {
	ROMol   *mol = (ROMol*)data;
	SparseFP	*res=NULL;
        std::vector<boost::uint32_t> invars(mol->getNumAtoms());
	try {
          RDKit::MorganFingerprints::getConnectivityInvariants(*mol,invars,true);
          res = (SparseFP*)RDKit::MorganFingerprints::getFingerprint(*mol, radius,&invars);
	} catch (...) {
		elog(ERROR, "makeMorganSFP: Unknown exception");
	}
	
	return (MolSparseFingerPrint)res;
}


extern "C" MolBitmapFingerPrint
makeMorganBFP(CROMol data, int radius) {
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;
        std::vector<boost::uint32_t> invars(mol->getNumAtoms());
	try {
          RDKit::MorganFingerprints::getConnectivityInvariants(*mol,invars,true);
          res = RDKit::MorganFingerprints::getFingerprintAsBitVect(*mol, radius,MORGAN_FP_SIZE,&invars);
	} catch (...) {
		elog(ERROR, "makeMorganBFP: Unknown exception");
	}
	
	return (MolBitmapFingerPrint)res;
}


extern "C" MolSparseFingerPrint 
makeAtomPairSFP(CROMol data){
	ROMol   *mol = (ROMol*)data;
	SparseFP	*res=NULL;
#ifdef UNHASHED_PAIR_FPS
	try {
          SparseIntVect<boost::int32_t> *afp=RDKit::AtomPairs::getAtomPairFingerprint(*mol);
          res = new SparseFP(1<<RDKit::AtomPairs::numAtomPairFingerprintBits);
          for(SparseIntVect<boost::int32_t>::StorageType::const_iterator iter=afp->getNonzeroElements().begin();
              iter!=afp->getNonzeroElements().end();++iter){
            res->setVal(iter->first,iter->second);
          }
          delete afp;
	} catch (...) {
		elog(ERROR, "makeAtomPairSFP: Unknown exception");
	}
#else
	try {
          SparseIntVect<boost::int32_t> *afp=RDKit::AtomPairs::getHashedAtomPairFingerprint(*mol,HASHED_PAIR_FP_SIZE);
          res = new SparseFP(HASHED_PAIR_FP_SIZE);
          for(SparseIntVect<boost::int32_t>::StorageType::const_iterator iter=afp->getNonzeroElements().begin();
              iter!=afp->getNonzeroElements().end();++iter){
            res->setVal(iter->first,iter->second);
          }
          delete afp;
	} catch (...) {
		elog(ERROR, "makeAtomPairSFP: Unknown exception");
	}
#endif	
	return (MolSparseFingerPrint)res;
}

extern "C" MolSparseFingerPrint 
makeTopologicalTorsionSFP(CROMol data){
	ROMol   *mol = (ROMol*)data;
	SparseFP	*res=NULL;

#ifdef UNHASHED_PAIR_FPS
	try {
          SparseIntVect<boost::int64_t> *afp=RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(*mol,boost::integer_traits<boost::uint32_t>::const_max);
          res = new SparseFP(boost::integer_traits<boost::uint32_t>::const_max);
          for(SparseIntVect<boost::int64_t>::StorageType::const_iterator iter=afp->getNonzeroElements().begin();
              iter!=afp->getNonzeroElements().end();++iter){
            res->setVal(iter->first,iter->second);
          }
          delete afp;
	} catch (...) {
		elog(ERROR, "makeTopologicalTorsionSFP: Unknown exception");
	}
#else
	try {
          SparseIntVect<boost::int64_t> *afp=RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(*mol,HASHED_PAIR_FP_SIZE);
          res = new SparseFP(HASHED_PAIR_FP_SIZE);
          for(SparseIntVect<boost::int64_t>::StorageType::const_iterator iter=afp->getNonzeroElements().begin();
              iter!=afp->getNonzeroElements().end();++iter){
            res->setVal(iter->first,iter->second);
          }
          delete afp;
	} catch (...) {
		elog(ERROR, "makeTopologicalTorsionSFP: Unknown exception");
	}
#endif
	return (MolSparseFingerPrint)res;
}

extern "C" MolBitmapFingerPrint 
makeAtomPairBFP(CROMol data){
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;
	try {
          res=RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(*mol,HASHED_PAIR_FP_SIZE);
	} catch (...) {
		elog(ERROR, "makeAtomPairBFP: Unknown exception");
	}
	return (MolBitmapFingerPrint)res;
}

extern "C" MolBitmapFingerPrint 
makeTopologicalTorsionBFP(CROMol data){
	ROMol   *mol = (ROMol*)data;
	ExplicitBitVect	*res=NULL;
	try {
          res =RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(*mol,HASHED_PAIR_FP_SIZE);
	} catch (...) {
          elog(ERROR, "makeTopologicalTorsionBFP: Unknown exception");
	}
	return (MolBitmapFingerPrint)res;
}

