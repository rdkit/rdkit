#include <string>

using namespace std;

#ifndef _Data3Ddescriptors

#define _Data3Ddescriptors

class Data3Ddescriptors
{
   private:
		static double mw[53];
		static double vdW[83];
		static double neg[53];
		static double pol[53];
		static double ionpol[83];
		static double rcov[83];

   public:
	   Data3Ddescriptors();
	   double* getMW();
	   double* getVDW();
	   double* getNEG();
	   double* getPOL();
	   double* getIonPOL();
	   double* getRCOV();

};
#endif 

