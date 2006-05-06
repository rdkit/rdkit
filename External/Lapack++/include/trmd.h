//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_TRIDIAG_MAT_DOUBLE_
#define _LA_TRIDIAG_MAT_DOUBLE_

#include "lafnames.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"

class LaTridiagMatDouble
{   
        LaVectorDouble du2_;    /* second upper diag, N-2 */
        LaVectorDouble du_;     /* upper diag, N-1 */
        LaVectorDouble d_;      /* main diag, N */
        LaVectorDouble dl_;     /* lower diag, N-1 */
        int size_;

        static double outofbounds_; /* return this address, when addresing out
                                    of bounds */
        static int debug_;        // print debug info.
        static int *info_;        // print matrix info only, not values
                                  //   originally 0, set to 1, and then
                                  //   reset to 0 after use.

    public:

        // constructors / destructor
    
        inline LaTridiagMatDouble();
        inline LaTridiagMatDouble(int N);
        inline LaTridiagMatDouble(const LaTridiagMatDouble &);
        inline ~LaTridiagMatDouble();


        // operators and member functions

            double & operator()(int i, int j);
            const double & operator()(int i, int j) const;
            LaVectorDouble diag(int); /* 0 main, -1 lower, 1 upper, 
                                            2 second upper  */
             LaVectorDouble diag(int) const; /* 0 main, -1 lower, 
                                            1 upper, 2 second upper  */
        inline LaTridiagMatDouble& ref(LaTridiagMatDouble&); 
        inline LaTridiagMatDouble& copy(const LaTridiagMatDouble&); 
        const LaTridiagMatDouble& info() const {
            int *t = info_; *t = 1; return *this;}
        int debug() const { return debug_;}
        int size() { return size_;}
        int size() const { return size_;}

        friend std::ostream& operator<<(std::ostream&,const LaTridiagMatDouble&);


};

    // constructors

inline LaTridiagMatDouble::LaTridiagMatDouble(): 
        du2_(), du_(), d_(), dl_(), size_(0)
{}

inline LaTridiagMatDouble::LaTridiagMatDouble(int N): 
        du2_(N-2), du_(N-1), d_(N), dl_(N-1), size_(N)
{}
    
inline LaTridiagMatDouble::LaTridiagMatDouble(const LaTridiagMatDouble& td):
        du2_(td.du2_), du_(td.du_), d_(td.d_), dl_(td.dl_), size_(td.size_)
{}

    // destructor

inline LaTridiagMatDouble::~LaTridiagMatDouble()
{
}


    // operators and member functions


inline LaTridiagMatDouble& LaTridiagMatDouble::ref(LaTridiagMatDouble&T) 
{
    du2_.ref(T.du2_);
    du_.ref(T.du_);
    d_.ref(T.d_);
    dl_.ref(T.dl_); 
    size_ = T.size_;

    return *this;
}


inline LaTridiagMatDouble& LaTridiagMatDouble::copy(const LaTridiagMatDouble&T) 
{
    du2_.copy(T.du2_);
    du_.copy(T.du_);
    d_.copy(T.d_);
    dl_.copy(T.dl_);    
    size_ = T.size_;

    return *this;
}





#endif 
// _LA_TRIDIAG_MAT_DOUBLE_
