//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.


#include "lafnames.h"
//#include LA_INDEX_H 
#include "laindex.h"

std::istream& operator>>(std::istream& s, LaIndex &i)
{
// NOTE: this only parses (s,e) indices.  Need to add code to 
//       or () case.
//
    char c;

    while ( (c = s.get()) != '(');      // ugly, but works...
    s >> i.start();                      // doesn't handle errors, though,
    while ( (c = s.get()) != ',');
    i.inc() = 1;
    s >> i.end();                        // no longer support inc() other
                                        // than inc == 1.
    while ( (c = s.get()) != ')'); 

    //  there's probably a way to do this with stream manipulators...
    //
    // s >> c >> i.start_ >> c >> i.inc_ >> c >> i.end_ >> c;
    //      (              ,              ,              ) `

    return s;
}
