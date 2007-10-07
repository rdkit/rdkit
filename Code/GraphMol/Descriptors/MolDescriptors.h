//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef _RD_MOLDESCRIPTORS_H_
#define _RD_MOLDESCRIPTORS_H_

namespace RDKit{
  class ROMol;
  namespace Descriptors {
    //! generate Wildman-Crippen LogP and MR estimates for a molecule
    /*!
      Uses an atom-based scheme based on the values in the paper:
        S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

      \param mol        the molecule of interest
      \param logp       used to return the logp estimate
      \param mr         used to return the MR estimate
      \param includeHs  (optional) if this is true (the default), a
          copy of \c mol is made and Hs are added to it.  If false,
	  Hs that are not explicitly present in the graph will not
	  be included.
	  
    */
    void CalcCrippenDescriptors(const ROMol &mol,double &logp,double &mr,
				bool includeHs=true);
    /*!
      Calculates a molecule's molecular weight

      \param mol        the molecule of interest
      \param onlyHeavy  (optional) if this is true (the default is false),
          only heavy atoms will be included in the MW calculation

      \return the AMW
    */
    double CalcAMW(const ROMol &mol,bool onlyHeavy=false);


  }
}

#endif
