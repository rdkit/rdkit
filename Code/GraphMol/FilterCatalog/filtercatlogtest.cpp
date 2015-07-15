//  $Id$
// 
//   Copyright (C) 2002-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>

#include <iostream>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

using namespace RDKit;
using namespace std;
RWMol _t;
typedef class ROMol Mol;

void testFilterCatalog() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2313979: aromaticity assignment hangs " << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    SmilesMolSupplier suppl(pathName+"pains.smi");

    FilterCatalogParams params;
    params.addCatalog(FilterCatalogParams::PAINS_A);
    params.addCatalog(FilterCatalogParams::PAINS_B);
    params.addCatalog(FilterCatalogParams::PAINS_C);
    
    FilterCatalog catalog(params);
    boost::scoped_ptr<ROMol> mol;
    const IntPair match1[] = {{0,23},{1,22},{2,20},{3,19},{4,25},{5,24},
                              {6,18},{7,17},{8,16},{9,21}};
    const IntPair match2[] = {{0,13},{1,12},{2,11},{3,14},{4,15},{5,10},
                              {6,9},{7,8},{8,7},{9,6},{10,5},{11,17},{12,16}};
    const IntPair match3[] = {{0,0},{1,1},{2,2},{3,4},{4,5},{5,6},{6,7},
                              {7,8},{8,9},{9,14},{10,15},{11,16}};
    int count = 0;
    while(!suppl.atEnd()){
      mol.reset(suppl.next());
      std::string name = mol->getProp<std::string>(common_properties::_Name);

      TEST_ASSERT(mol.get());
      if (catalog.hasMatch(*mol)) {
        std::cerr << "Warning: molecule failed filter " << std::endl;
      }
      // More detailed
      FilterCatalog::CONST_SENTRY entry = catalog.getFirstMatch(*mol);
      if (entry) {
        std::cerr << "Warning: molecule failed filter: reason " <<
          entry->getDescription() << std::endl;
        TEST_ASSERT(entry->getDescription() == name);
        
        // get the substructure atoms for visualization
        std::vector<FilterMatch> matches;
        if (entry->getFilterMatches(*mol, matches)) {
          for(std::vector<FilterMatch>::const_iterator it = matches.begin();
              it != matches.end(); ++it) {
            // Get the FilterMatcherBase that matched
            const FilterMatch & fm = (*it);
            boost::shared_ptr<FilterMatcherBase> matchingFilter =       \
              fm.filterMatch;
            
            // Get the matching atom indices
            const MatchVectType &vect = fm.atomPairs;
            switch(count) {
            case 0: TEST_ASSERT(check(vect, &match1[0])); break;
            case 1: TEST_ASSERT(check(vect, &match2[0])); break;
            case 2: TEST_ASSERT(check(vect, &match3[0])); break;
            }
              
            // do something with these...
          }
        }
      }
      count ++;
    } // end while
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
  

int main(){
  RDLog::InitLogs();
  //boost::logging::enable_logs("rdApp.debug");


  testFilterCatalog();

  return 0;
}


