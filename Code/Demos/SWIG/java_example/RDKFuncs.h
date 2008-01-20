#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>

RDKit::ROMol *MolFromSmiles(std::string smi);
RDKit::ROMol *MolFromSmarts(std::string sma);
std::string MolToSmiles(RDKit::ROMol *mol,bool doIsomericSmiles=false,
                        bool doKekule=false, int rootedAtAtom=-1);
bool HasSubstructMatch(RDKit::ROMol &mol,RDKit::ROMol &query,bool useChirality=false,
                       bool registerQuery=false);

