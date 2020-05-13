#include <RDGeneral/test.h>
#include <cstdlib>
#include <cstdio>

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>

#define VERBOSE

static void testToPDB(const RDKit::ROMol &mol) {
  std::string pdb = RDKit::MolToPDBBlock(mol);
  printf("  PDB:\n%s\n", pdb.c_str());
}

static void testToSmiles(const RDKit::ROMol &mol) {
  std::string smi = RDKit::MolToSmiles(mol, true);
  printf("  Smiles: %s\n", smi.c_str());
}

static void testToSequence(const RDKit::ROMol &mol) {
  std::string seq = RDKit::MolToSequence(mol);
  printf("  Sequence: %s\n", seq.c_str());
}

static void testToHELM(const RDKit::ROMol &mol) {
  std::string helm = RDKit::MolToHELM(mol);
  printf("  HELM: %s\n", helm.c_str());
}

static void testSeq(const char *seq) {
  RDKit::RWMol *mol = RDKit::SequenceToMol(seq);
  testToSequence(*mol);
  testToHELM(*mol);
  testToSmiles(*mol);
  testToPDB(*mol);
  delete mol;
}

static void testHELM(const char *helm) {
  RDKit::RWMol *mol = RDKit::HELMToMol(helm);
  TEST_ASSERT(mol);
  testToSequence(*mol);
  testToHELM(*mol);
  testToSmiles(*mol);
  testToPDB(*mol);
  delete mol;
}

static void testPDB(const char *pdb) {
  RDKit::RWMol *mol = RDKit::PDBBlockToMol(pdb);
  testToSequence(*mol);
  testToHELM(*mol);
  testToSmiles(*mol);
  testToPDB(*mol);
  delete mol;
}

static void testPDB() {
  testPDB(
      "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      3  C   GLY A   1       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      4  O   GLY A   1       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      5  N   ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      6  CA  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      7  C   ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      8  O   ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM      9  CB  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     10  CG  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     11  CD  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     12  NE  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     13  CZ  ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     14  NH1 ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     15  NH2 ARG A   2       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     16  N   GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     17  CA  GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     18  C   GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     19  O   GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     20  CB  GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     21  CG  GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     22  CD  GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     23  OD1 GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     24  OD2 GLU A   3       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     25  N   GLY A   4       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     26  CA  GLY A   4       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     27  C   GLY A   4       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     28  O   GLY A   4       0.000   0.000   0.000  1.00  0.00\n"
      "ATOM     29  OXT GLY A   4       0.000   0.000   0.000  1.00  0.00\n"
      "CONECT    1    2\n"
      "CONECT    2    1    3\n"
      "CONECT    3    2    4    4    5\n"
      "CONECT    4    3    3\n"
      "CONECT    5    3    6\n"
      "CONECT    6    5    7    9\n"
      "CONECT    7    6    8    8   16\n"
      "CONECT    8    7    7\n"
      "CONECT    9    6   10\n"
      "CONECT   10    9   11\n"
      "CONECT   11   10   12\n"
      "CONECT   12   11   13\n"
      "CONECT   13   12   14   14   15\n"
      "CONECT   14   13   13\n"
      "CONECT   15   13\n"
      "CONECT   16    7   17\n"
      "CONECT   17   16   18   20\n"
      "CONECT   18   17   19   19   25\n"
      "CONECT   19   18   18\n"
      "CONECT   20   17   21\n"
      "CONECT   21   20   22\n"
      "CONECT   22   21   23   23   24\n"
      "CONECT   23   22   22\n"
      "CONECT   24   22\n"
      "CONECT   25   18   26\n"
      "CONECT   26   25   27\n"
      "CONECT   27   26   28   28   29\n"
      "CONECT   28   27   27\n"
      "CONECT   29   27\n"
      "END\n");
}

int main() {
  testSeq("GREG");

  // The following is the sequence of crambin from 1CRN
  testSeq("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN");

  testHELM("PEPTIDE1{G.R.E.G}$$$$");

  // The following is the HELM for oxytocin
  testHELM("PEPTIDE1{C.Y.I.Q.N.C.P.L.G.[am]}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$");

  // The following is the HELM for gramicidin S
  testHELM(
      "PEPTIDE1{L.[dF].P.V.[Orn].L.[dF].P.V.[Orn]}$"
      "PEPTIDE1,PEPTIDE1,10:R2-1:R1$$$");

  // The following tests GREG from PDB file
  testPDB();

  return 0;
}
