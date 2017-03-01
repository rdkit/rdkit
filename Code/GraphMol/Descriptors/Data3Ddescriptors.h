#include <string>

using namespace std;

#ifndef _Data3Ddescriptors

#define _Data3Ddescriptors

class Data3Ddescriptors
{
   private:
		static double mw[110];
		static double vdW[110];
		static double neg[110];
		static double pol[110];
		static double ionpol[110];
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

