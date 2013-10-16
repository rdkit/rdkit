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
    fputs("usage:  smi2pdb <infile.smi> [<outfile.pdb>]\n",stderr);
    exit(1);
  }

  // filenName, delimiter=" \t", smilesColumn=0, nameColumn=1
  // titleLine=true, sanitize=true
  RDKit::SmilesMolSupplier supplier(inpname," \t",0,1,false,true);
  RDKit::PDBWriter writer(outname);
  while (!supplier.atEnd()) {
    RDKit::ROMol *mol = supplier.next();
    if (mol) {
      writer.write(*mol);
      delete mol;
    }
  }
  writer.close();
  return 0;
}

