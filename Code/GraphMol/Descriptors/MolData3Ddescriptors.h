#include <string>
#include "Data3Ddescriptors.h"
#include <GraphMol/RDKitBase.h>


using namespace std;

#ifndef _MolData3Ddescriptors

#define _MolData3Ddescriptors

class MolData3Ddescriptors
{
	private:
		Data3Ddescriptors data3D;


	public:
	   MolData3Ddescriptors();
	   std::vector<double> GetCharges(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativeMW(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativePol(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativeRcov(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativeENeg(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativeIonPol(const RDKit::ROMol& mol);
	   std::vector<double> GetRelativeVdW(const RDKit::ROMol& mol);
	   std::vector<double> GetUn(int numAtoms);
	   int GetPrincipalQuantumNumber(int AtomicNum);
	   std::vector<double> GetIState(const  RDKit::ROMol& mol);
	   std::vector<double> GetEState(const  RDKit::ROMol& mol);
	   std::vector<double> GetEState2(const  RDKit::ROMol& mol);


};
#endif 


