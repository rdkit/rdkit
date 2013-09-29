#include <stdlib.h>
#include <stdio.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

int main(int argc, char *argv[])
{
  const char *inpname;
  const char *outname;

  if (argc == 2) {
    inpname = argv[1];
    outname = "-";
  } else if (argc == 3) {
    inpname = argv[1];
    outname = argv[2];
  } else {
    fputs("usage:  pdb2smi <infile.pdb> [<outfile.smi>]\n",stderr);
    exit(1);
  }

  // fileName, sanitize=true, removeHs=true, strictParsing=true
  RDKit::PDBMolSupplier supplier(inpname,true,true,0);
  // fileName, delimiter=" ", nameHeader="Name"
  // includeHeader=true, isomericSmiles=false, kekuleSmiles=false
  RDKit::SmilesWriter writer(outname," ","Name",false,true,false);
  while (!supplier.atEnd()) {
    RDKit::ROMol *mol = supplier.next();
    if (mol) {
      writer.write(*mol);
      delete mol;
    } else std::cout << "#" << std::endl;
  }
  writer.close();
  return 0;
}

