// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//       All Rights Reserved
//
// The tests in this file demonstrate the use of RDKit with
// the boost iostreams library. Particularly the compression pieces
// of that library.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <fstream>
#include <iostream>
#ifdef SUPPORT_COMPRESSED_IO
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/range/iterator_range.hpp>
#endif
#include <string>

using namespace RDKit;
namespace io=boost::iostreams;

void test1(){
#ifdef SUPPORT_COMPRESSED_IO
  BOOST_LOG(rdInfoLog) << "testing basic reading from compressed streams" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fname2 = rdbase + "/Code/Demos/RDKit/BinaryIO/test_data/triazine.mol.gz";
  RWMol *m;
  
  io::filtering_istream inStrm;
  inStrm.push(io::gzip_decompressor());
  inStrm.push(io::file_source(fname2));
  TEST_ASSERT(inStrm.is_complete());
  unsigned int lineNo=0;
  m = MolDataStreamToMol(inStrm,lineNo);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(6));

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
#else
  BOOST_LOG(rdInfoLog) << "Compressed IO disabled; test skipped" << std::endl;
#endif
}

void test2(){
#ifdef SUPPORT_COMPRESSED_IO
  BOOST_LOG(rdInfoLog) << "testing writing to, then reading from compressed streams" << std::endl;

  std::string smiles="C1CCCC1";
  std::string buff,molBlock;
  RWMol *m;
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(5));
  m->setProp("_Name","monkey");
  io::filtering_ostream outStrm;
  outStrm.push(io::gzip_compressor());
  outStrm.push(io::back_inserter(buff));
  TEST_ASSERT(outStrm.is_complete());

  molBlock = MolToMolBlock(*m);
  outStrm << molBlock;
  outStrm.reset();
  
  delete m;
  io::filtering_istream inStrm;
  inStrm.push(io::gzip_decompressor());
  inStrm.push(boost::make_iterator_range(buff));
  TEST_ASSERT(inStrm.is_complete());

  unsigned int lineNo=0;
  m = MolDataStreamToMol(inStrm,lineNo);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(5));

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
#else
  BOOST_LOG(rdInfoLog) << "Compressed IO disabled; test skipped" << std::endl;
#endif
}


void test3(){
#ifdef SUPPORT_COMPRESSED_IO
  BOOST_LOG(rdInfoLog) << "testing writing pickles to a file then reading them back" << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fname2 = rdbase + "/Code/Demos/RDKit/BinaryIO/test_data/mols.rdb";

  std::string smiles,buff;
  RWMol *m;

  std::vector<unsigned int> filePs;
  io::filtering_ostream outStrm;
  outStrm.push(io::file_sink(fname2,std::ios_base::out|std::ios_base::binary));
  TEST_ASSERT(outStrm.is_complete());
  
  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));
  
  MolPickler::pickleMol(*m,outStrm);
  delete m;

  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(0);
  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="C1CCCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(5));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="c1ccccc1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(6));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="c1ccccc1CC(=O)O";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(9));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,outStrm);
  delete m;

  outStrm.flush();
  outStrm.reset();
  
  io::filtering_istream inStrm;
  inStrm.push(io::file_source(fname2,std::ios_base::in|std::ios_base::binary));
  TEST_ASSERT(inStrm.is_complete());

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(5));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(6));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(9));
  delete m;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
#else
  BOOST_LOG(rdInfoLog) << "Compressed IO disabled; test skipped" << std::endl;
#endif
}

void test4(){
#ifdef SUPPORT_COMPRESSED_IO
  BOOST_LOG(rdInfoLog) << "testing writing pickles to a single compressed file then reading them back" << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fname2 = rdbase + "/Code/Demos/RDKit/BinaryIO/test_data/mols.rdz";

  std::string smiles;
  RWMol *m;

  io::filtering_ostream outStrm;
  outStrm.push(io::gzip_compressor());
  outStrm.push(io::file_sink(fname2,std::ios_base::out|std::ios_base::binary));
  TEST_ASSERT(outStrm.is_complete());
  
  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));
  
  MolPickler::pickleMol(*m,outStrm);
  delete m;

  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));
  RDDepict::compute2DCoords(*m);

  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="C1CCCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(5));
  RDDepict::compute2DCoords(*m);
  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="c1ccccc1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(6));
  RDDepict::compute2DCoords(*m);
  MolPickler::pickleMol(*m,outStrm);
  delete m;
  
  smiles="c1ccccc1CC(=O)O";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(9));
  RDDepict::compute2DCoords(*m);
  MolPickler::pickleMol(*m,outStrm);
  delete m;

  io::flush(outStrm);
  outStrm.pop();

  
  io::filtering_istream inStrm;
  inStrm.push(io::gzip_decompressor());
  inStrm.push(io::file_source(fname2,std::ios_base::in|std::ios_base::binary));
  TEST_ASSERT(inStrm.is_complete());

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(5));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(6));
  delete m;

  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(inStrm,*m);
  TEST_ASSERT(m->getNumAtoms(9));
  delete m;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
#else
  BOOST_LOG(rdInfoLog) << "Compressed IO disabled; test skipped" << std::endl;
#endif
}

void test5(){
#ifdef SUPPORT_COMPRESSED_IO
  BOOST_LOG(rdInfoLog) << "testing writing compressed pickles to a single file then reading them back" << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fname2 = rdbase + "/Code/Demos/RDKit/BinaryIO/test_data/tmp.rdz";

  std::string smiles,buff;
  RWMol *m;
  std::vector<unsigned int> filePs;

  io::filtering_ostream outStrm;
  outStrm.push(io::file_sink(fname2,std::ios_base::out|std::ios_base::binary));
  TEST_ASSERT(outStrm.is_complete());
  
  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));


  io::filtering_ostream *tmpStrm;
  tmpStrm=new io::filtering_ostream();
  tmpStrm->push(io::gzip_compressor());
  tmpStrm->push(io::back_inserter(buff));
  
  
  filePs.push_back(0);
  MolPickler::pickleMol(*m,*tmpStrm);
  delete m;
  tmpStrm->reset();

  outStrm<<buff.size();
  outStrm<<buff;
  std::cerr<<"sz: " <<buff.size()<<" "<<outStrm.tellp() <<std::endl;
  buff="";
  tmpStrm->push(io::gzip_compressor());
  tmpStrm->push(io::back_inserter(buff));
  
  smiles="C1CCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(4));
  RDDepict::compute2DCoords(*m);

  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,*tmpStrm);
  delete m;
  tmpStrm->reset();
  outStrm<<buff.size();
  outStrm<<buff;
  std::cerr<<"sz: " <<buff.size()<<" "<<outStrm.tellp() <<std::endl;
  buff="";
  tmpStrm->push(io::gzip_compressor());
  tmpStrm->push(io::back_inserter(buff));

  smiles="C1CCCC1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(5));
  RDDepict::compute2DCoords(*m);

  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,*tmpStrm);
  delete m;
  tmpStrm->reset();
  outStrm<<buff.size();
  outStrm<<buff;
  std::cerr<<"sz: " <<buff.size()<<" "<<outStrm.tellp() <<std::endl;
  buff="";
  tmpStrm->push(io::gzip_compressor());
  tmpStrm->push(io::back_inserter(buff));


  smiles="c1ccccc1";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(6));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,*tmpStrm);
  delete m;
  tmpStrm->reset();
  outStrm<<buff.size();
  outStrm<<buff;
  std::cerr<<"sz: " <<buff.size()<<" "<<outStrm.tellp() <<std::endl;
  buff="";
  tmpStrm->push(io::gzip_compressor());
  tmpStrm->push(io::back_inserter(buff));
  
  smiles="c1ccccc1CC(=O)O";
  m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms(9));
  RDDepict::compute2DCoords(*m);
  filePs.push_back(outStrm.tellp());
  MolPickler::pickleMol(*m,*tmpStrm);
  delete m;
  tmpStrm->reset();
  outStrm<<buff.size();
  outStrm<<buff;
  std::cerr<<"sz: " <<buff.size()<<" "<<outStrm.tellp() <<std::endl;

  delete tmpStrm;
  io::flush(outStrm);
  outStrm.pop();

  io::filtering_istream tmpIStrm;
  io::filtering_istream inStrm;
  unsigned int sz;
  char *charArr;
  
  inStrm.push(io::file_source(fname2,std::ios_base::in|std::ios_base::binary));
  TEST_ASSERT(inStrm.is_complete());

  inStrm>>sz;
  charArr = new char[sz];
  inStrm.read(charArr,sz);
  buff = "";
  buff.append(charArr,sz);
  tmpIStrm.push(io::gzip_decompressor());
  tmpIStrm.push(boost::make_iterator_range(buff));
  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(tmpIStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  inStrm>>sz;
  delete [] charArr;
  charArr = new char[sz];
  inStrm.read(charArr,sz);
  buff = "";
  buff.append(charArr,sz);
  tmpIStrm.reset();
  tmpIStrm.push(io::gzip_decompressor());
  tmpIStrm.push(boost::make_iterator_range(buff));
  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(tmpIStrm,*m);
  TEST_ASSERT(m->getNumAtoms(4));
  delete m;

  inStrm>>sz;
  delete [] charArr;
  charArr = new char[sz];
  inStrm.read(charArr,sz);
  buff = "";
  buff.append(charArr,sz);
  tmpIStrm.reset();
  tmpIStrm.push(io::gzip_decompressor());
  tmpIStrm.push(boost::make_iterator_range(buff));
  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(tmpIStrm,*m);
  TEST_ASSERT(m->getNumAtoms(5));
  delete m;

  inStrm>>sz;
  delete [] charArr;
  charArr = new char[sz];
  inStrm.read(charArr,sz);
  buff = "";
  buff.append(charArr,sz);
  tmpIStrm.reset();
  tmpIStrm.push(io::gzip_decompressor());
  tmpIStrm.push(boost::make_iterator_range(buff));
  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(tmpIStrm,*m);
  TEST_ASSERT(m->getNumAtoms(6));
  delete m;

  inStrm>>sz;
  delete [] charArr;
  charArr = new char[sz];
  inStrm.read(charArr,sz);
  buff = "";
  buff.append(charArr,sz);
  tmpIStrm.reset();
  tmpIStrm.push(io::gzip_decompressor());
  tmpIStrm.push(boost::make_iterator_range(buff));
  m = new RWMol();
  TEST_ASSERT(m);
  MolPickler::molFromPickle(tmpIStrm,*m);
  TEST_ASSERT(m->getNumAtoms(9));
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
#else
  BOOST_LOG(rdInfoLog) << "Compressed IO disabled; test skipped" << std::endl;
#endif
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
#endif
  return 0;
}
