//
// Copyright (C)  2001-2013 Greg Landrum, Randal M. Henne and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef __RD_INVARIANT_H__
#define __RD_INVARIANT_H__

#include <assert.h>
#include <string>
#include <iostream>
#include <stdexcept>

#include <RDGeneral/RDLog.h>
// 
// What if no invariant method is defined?
//
#if !defined INVARIANT_EXCEPTION_METHOD && \
    !defined INVARIANT_ASSERT_METHOD && \
    !defined INVARIANT_SILENT_METHOD
#define INVARIANT_EXCEPTION_METHOD 1
#endif

//
// What if an invariant method is defined, but none are true?
//
#if !INVARIANT_EXCEPTION_METHOD && \
    !INVARIANT_ASSERT_METHOD && \
    !INVARIANT_SILENT_METHOD
#undef  INVARIANT_EXCEPTION_METHOD
#define INVARIANT_EXCEPTION_METHOD 1
#endif


namespace Invar {

  class Invariant : public std::runtime_error {

  public:
  
    Invariant( const char * prefix, const char * mess, const char * expr, const char * const file, int line )
      : std::runtime_error( prefix ),
	mess_d( mess ),
        expr_d( expr ),
	prefix_d( prefix ),
	file_dp( file ),
	line_d( line )
    {
    }
    Invariant( const char * prefix, std::string mess, const char * expr, const char * const file, int line )
      : std::runtime_error( prefix ),
	mess_d( mess.c_str() ),
        expr_d( expr ),
	prefix_d( prefix ),
	file_dp( file ),
	line_d( line )
    {
    }
    ~Invariant() throw () {};

    std::string 
    getMessage() const
    { return mess_d; }

    const char *
    getFile() const
    { return file_dp; }

    std::string
    getExpression() const
    { return expr_d; }

    int
    getLine() const
    { return line_d; }

    std::string  
    toString() const;
    

  private:

    std::string
      mess_d,
      expr_d,prefix_d;
    

    
    const char 
      * const file_dp;

    int
      line_d;
  };
  std::ostream& operator<<( std::ostream & s, const Invariant & inv );
} // end of namespace Invar



#define ASSERT_INVARIANT( expr, mess ) \
    assert( expr )

//
// Set desired reporting method
//

#if INVARIANT_EXCEPTION_METHOD

#define CHECK_INVARIANT( expr, mess ) if ( !(expr) ) {\
     Invar::Invariant inv( "Invariant Violation", mess, \
     #expr, __FILE__, __LINE__ ); \
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv; }

#define PRECONDITION( expr, mess ) if ( !(expr) ) {\
     Invar::Invariant inv( "Pre-condition Violation", mess, \
     #expr, __FILE__, __LINE__ ); \
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv; }

#define POSTCONDITION( expr, mess ) if ( !(expr) ) { \
     Invar::Invariant inv( "Post-condition Violation", mess, \
     #expr, __FILE__, __LINE__ );\
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv; }

#define UNDER_CONSTRUCTION( fn ) Invar::Invariant inv( "Incomplete Code", \
     "This routine is still under development", fn, __FILE__, __LINE__ ); \
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv;

#define RANGE_CHECK( lo, x, hi ) if ( (lo)>(hi) || (x)<(lo) || (x)>(hi) ) {\
     std::stringstream errstr;\
     errstr << lo << " <= " << x <<" <= "<<hi;\
     Invar::Invariant inv( "Range Error", #x, errstr.str().c_str(), __FILE__, __LINE__ );\
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv; }

#define TEST_ASSERT( expr ) if ( !(expr) ) {\
     Invar::Invariant inv( "Test Assert", "Expression Failed: ", \
     #expr, __FILE__, __LINE__ ); \
     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n\n"; throw inv; }


#elif INVARIANT_ASSERT_METHOD

#define CHECK_INVARIANT( expr, mess ) assert( expr );
#define PRECONDITION( expr, mess ) assert( expr );
#define POSTCONDITION( expr, mess ) assert( expr );
#define UNDER_CONSTRUCTION(fn ) assert(0);
#define RANGE_CHECK( lo, x, hi ) assert( (lo)<=(hi) && (x)>=(lo) && (x)<=(hi) );
#define TEST_ASSERT( expr ) assert(expr);

#elif INVARIANT_SILENT_METHOD

#define CHECK_INVARIANT( expr, mess ) 
#define PRECONDITION( expr, mess ) 
#define POSTCONDITION( expr, mess )
#define UNDER_CONSTRUCTION( fn )
#define RANGE_CHECK( lo, x, hi )
#define TEST_ASSERT( expr )

#endif

#endif

