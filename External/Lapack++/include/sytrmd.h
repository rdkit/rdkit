//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
#define _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_

//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"

class LaSymmTridiagMatDouble
{   
        int size_;
        LaVectorDouble d_;          /* main diag */
        LaVectorDouble dl_;         /* lower diag */
        static double outofbounds_; /* return this address, when addresing 
                                        out of bounds */
        static int debug_;        // print debug info.
        static int *info_;        // print matrix info only, not values
                                  //   originally 0, set to 1, and then
                                  //   reset to 0 after use.

    public:

        // constructors / destructor
    
        inline LaSymmTridiagMatDouble();
        inline LaSymmTridiagMatDouble(int N);
        inline LaSymmTridiagMatDouble(const LaSymmTridiagMatDouble &);
        inline ~LaSymmTridiagMatDouble();

        // operators and member functions

                double& operator()(int i, int j);
                double operator()(int i, int j) const;
                LaVectorDouble diag(int); /* 0 main, -1 lower, 1 upper */
                LaVectorDouble diag(int) const;
        inline LaSymmTridiagMatDouble& ref(LaSymmTridiagMatDouble&); 
        inline LaSymmTridiagMatDouble& copy(const LaSymmTridiagMatDouble&); 
        inline const LaSymmTridiagMatDouble& info() const {
            int *t = info_;
            *t = 1;
            return *this;};
        inline int debug() const {    // return debug flag.
            return debug_;}
        inline int size() const; 



        friend std::ostream& operator<<(std::ostream&,const LaSymmTridiagMatDouble&);


};

    // constructors

inline LaSymmTridiagMatDouble::LaSymmTridiagMatDouble(): 
        d_(), dl_()
{
    size_ = 0;
}

inline LaSymmTridiagMatDouble::LaSymmTridiagMatDouble(int N): 
        d_(N), dl_(N-1)
{
    size_ = N;
}
    
inline LaSymmTridiagMatDouble::LaSymmTridiagMatDouble(const LaSymmTridiagMatDouble& td): d_(td.d_), dl_(td.dl_)
{
    size_ = td.size_;
}

    // destructor

inline LaSymmTridiagMatDouble::~LaSymmTridiagMatDouble()
{
}


    // operators and member functions



inline LaSymmTridiagMatDouble& LaSymmTridiagMatDouble::ref(LaSymmTridiagMatDouble&T) 
{
    d_.ref(T.d_);
    dl_.ref(T.dl_); 
    size_ = T.size_;

    return *this;
}


inline LaSymmTridiagMatDouble& LaSymmTridiagMatDouble::copy(const LaSymmTridiagMatDouble&T) 
{
    d_.copy(T.d_);
    dl_.copy(T.dl_);    
    size_ = T.size_;

    return *this;
}

inline int LaSymmTridiagMatDouble::size() const
{
    return size_;
}



#endif 
// _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
