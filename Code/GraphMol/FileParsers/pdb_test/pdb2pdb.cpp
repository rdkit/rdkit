#include <stdlib.h>
#include <stdio.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/PDBSupplier.h>


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
    fputs("usage:  pdb2pdb <infile.pdb> [<outfile.pdb>]\n",stderr);
    exit(1);
  }

  // fileName, sanitize=true, removeHs=true, flavor=0
  RDKit::PDBMolSupplier supplier(inpname,true,true,0);
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

