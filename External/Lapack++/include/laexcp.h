//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

//
//  C++ exception handling is not currently supported by most compilers.
//  The macros below allow one to "mimic" the throw expression of
//  character strings.  Transfer of control does not correspond
//  to a "try" block, but rather exits the program with the proper
//  banner.  This is a similar behavior to the real exception handling
//  if there was no explicit user-supplied try block.
//

#ifndef _LA_EXCEPTION_H_
#define _LA_EXCEPTION_H_

#include <iostream>
#include <stdlib.h>

#define LaException(where, what)    where , what
#define throw throw_
inline void throw(const char *where, const char *what)
{
    std::cerr << "Exception: " << where << "  " << what << std::endl;
    exit(1);
}

#endif  

