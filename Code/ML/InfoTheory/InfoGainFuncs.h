// $Id$
//
//  Copyright (C) 2003 Rational Discovery LLC
//

#ifndef INFOGAINFUNC_H
#define INFOGAINFUNC_H

#include <RDGeneral/types.h>

namespace RDInfoTheory {

  template<class T> double ChiSquare(T *dMat, long int dim1,long int dim2) {
    // For a contingency matrix with each column corresponding to a class and each row to a 
    // the descriptor (or variable) state, the matrix looks something like for 3x3 problem
    // 
    //            1    2    3   Totals
    //      1 |  N11  N12  N13    R1
    //      2 |  N21  N22  N23    R2
    //      3 |  N31  N32  N33    R3
    // Totals |   C1   C2   C3    N
    //
    //  Th chi squere formula is 
    //  chi = sum((N/Ri)*sum(Nij^2/Cj) ) -N
    T *rowSums, *colSums;
    int i, j, tSum;
    // find the row sum
    tSum = 0;
    rowSums = new T[dim1];
    for (i = 0; i < dim1; i++) {
      int idx1 = i*dim2;
      rowSums[i] = (T)0.0;
      for (j = 0; j < dim2; j++) {
        rowSums[i] += dMat[idx1 + j];
      }
      tSum += (int)rowSums[i];
    }

    // find the column sums
    colSums = new T[dim2];
    for (i = 0; i < dim2; i++) {
      colSums[i] = (T)0.0;
      for (j = 0; j < dim1; j++) {
        colSums[i] += dMat[j*dim2 + i];
      }
    }
    
    double chi = 0.0;
    for ( i = 0; i < dim1; i++) {
      double rchi = 0.0;
      for (j = 0; j < dim2; j++) {
        rchi += (pow((double)dMat[i*dim2 + j], 2)/colSums[j]);
      }
      chi += ( ((double)tSum/rowSums[i])*rchi );
    }
    chi -= tSum;
    delete [] rowSums;
    delete [] colSums;

    return chi;
  }

  template<class T> double InfoEntropy(T *tPtr, long int dim) {
    int i;
    T nInstances = 0;
    double accum=0.0,d;
    
    for(i=0;i<dim;i++){
      nInstances += tPtr[i];
    }
  
    if(nInstances != 0){
      for(i=0;i<dim;i++){
        d = (double)tPtr[i]/nInstances;
        if(d != 0){
          accum += -d*log(d);
        }
      }
    }
    return accum/log(2.0);
  }

  template<class T> double InfoEntropyGain(T *dMat, long int dim1,long int dim2) {
    T *variableRes, *overallRes;
    double gain,term2;
    int tSum;

    //std::cerr<<" --------\n    ieg: "<<dim1<<" "<<dim2<<std::endl;
    variableRes = new T[dim1];
    for(long int i=0;i<dim1;i++){
      long int idx1 = i*dim2;
      variableRes[i] = (T)0.0;
      for(long int j=0;j<dim2;j++){
        variableRes[i] += dMat[idx1+j];
        //std::cerr<<"  "<<i<<" "<<j<<" "<<dMat[idx1+j]<<std::endl;
      }
    }

    overallRes = new T[dim2];
    // do the col sums
    for(long int i=0;i<dim2;i++){
      overallRes[i] = (T)0.0;
      for(long int j=0;j<dim1;j++){
        overallRes[i] += dMat[j*dim2+i];
        //std::cerr<<"  "<<i<<" "<<j<<" "<<dMat[j*dim2+i]<<std::endl;
      }
    }

    term2 = 0.0;
    for(long int i=0;i<dim1;i++) {
      T *tPtr;
      tPtr = dMat + i*dim2;
      term2 += variableRes[i] * InfoEntropy(tPtr,dim2);
    }
    tSum = 0;
    for(long int i=0;i<dim2;i++){
      tSum += static_cast<int>(overallRes[i]);
    }
    
    if(tSum != 0){
      term2 /= tSum;
      gain = InfoEntropy(overallRes,dim2) - term2;
    }
    else{
      gain = 0.0;
    }
    //std::cerr<<"  >gain> "<<gain<<std::endl;
    
    delete [] overallRes;
    delete [] variableRes;
    return gain;
  }
   
  
}
#endif


