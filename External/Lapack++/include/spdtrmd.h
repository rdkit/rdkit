//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SPD_TRIDIAG_MAT_DOUBLE_H_
#define _LA_SPD_TRIDIAG_MAT_DOUBLE_H_

#include "lafnames.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"

class LaSpdTridiagMatDouble
{   
        int size_;
        LaVectorDouble d_;  /* main diag */
        LaVectorDouble dl_; /* lower diag */
        static double outofbounds_; /* return this address, when addresing 
                                        out of bounds */
        static int debug_;        // print debug info.
        static int *info_;        // print matrix info only, not values
                                  //   originally 0, set to 1, and then
                                  //   reset to 0 after use.

    public:

        // constructors / destructor
    
        inline LaSpdTridiagMatDouble();
        inline LaSpdTridiagMatDouble(int N);
        inline LaSpdTridiagMatDouble(const LaSpdTridiagMatDouble &);
        inline ~LaSpdTridiagMatDouble();

        // operators and member functions

            double& operator()(int i, int j);
            double operator()(int i, int j) const;
            LaVectorDouble diag(int); /* 0 main, -1 lower, 1 upper */
        inline LaSpdTridiagMatDouble& ref(LaSpdTridiagMatDouble&); 
        inline LaSpdTridiagMatDouble& copy(const LaSpdTridiagMatDouble&); 
        inline const LaSpdTridiagMatDouble& info() const {
            int *t = info_;
            *t = 1;
            return *this;};
        inline int debug() const {    // return debug flag.
            return debug_;}
        inline int size() const; 



        friend std::ostream& operator<<(std::ostream&,const LaSpdTridiagMatDouble&);


};

    // constructors

inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble(): 
        d_(), dl_()
{
    size_ = 0;
}

inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble(int N): 
        d_(N), dl_(N-1)
{
    size_ = N;
}
    
inline LaSpdTridiagMatDouble::LaSpdTridiagMatDouble(const LaSpdTridiagMatDouble& td): d_(td.d_), dl_(td.dl_)
{
    size_ = td.size_;
}

    // destructor

inline LaSpdTridiagMatDouble::~LaSpdTridiagMatDouble()
{
}


    // operators and member functions





inline LaSpdTridiagMatDouble& LaSpdTridiagMatDouble::ref(LaSpdTridiagMatDouble&T) 
{
    d_.ref(T.d_);
    dl_.ref(T.dl_); 
    size_ = T.size_;

    return *this;
}


inline LaSpdTridiagMatDouble& LaSpdTridiagMatDouble::copy(const LaSpdTridiagMatDouble&T) 
{
    d_.copy(T.d_);
    dl_.copy(T.dl_);    
    size_ = T.size_;

    return *this;
}

inline int LaSpdTridiagMatDouble::size() const
{
    return size_;
}



#endif 
// _LA_SPD_TRIDIAG_MAT_DOUBLE_H_
