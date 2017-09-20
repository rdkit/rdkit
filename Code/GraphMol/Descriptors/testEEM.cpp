//  Created by Guillaume GODIN
//  "Copyright 2013-2016 Tomas Racek (tom@krab1k.net)"
//  Copyright (C) 2012-2017 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <GraphMol/Descriptors/EEM.h>

void testEEM() {
  std::cout << "=>start test EEM\n";
  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/set00.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string line;
  int nDone = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    int confId=-1;
    std::vector<double> charges;

    RDKit::Descriptors::EEM(*m, charges);

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testEEM();
}
