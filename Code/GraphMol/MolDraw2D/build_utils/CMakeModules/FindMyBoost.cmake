# Attempts to fine Boost libraries using the current value of $BOOST,
# set by module load boost on our system. It won't look in the default
# places, and it needs boost 1.40 or later (without the library name decorations
# that boost used to use).
# It also returns the static (.a) libraries not the .so ones because that's
# easiest for shipping (on Unix anyway. This may have to change once I start
# dealing with Vista as well.)
#
# It will define
# BOOST_FOUND as MyBoost_FOUND if it finds everything it needs
# BOOST_INCLUDE_DIR
# BOOST_LIBRARIES

set(BOOST_DIR $ENV{BOOST})
if(NOT BOOST_DIR)
  message( FATAL_ERROR "Environment variable BOOST must be set." )
endif(NOT BOOST_DIR)

# include dir
find_path(BOOST_INCLUDE_DIR
  boost/lexical_cast.hpp
  ${BOOST_DIR}/include NO_DEFAULT_PATH)

if(NOT BOOST_INCLUDE_DIR)
  message( FATAL_ERROR "Couldn't find boost include directory.")
endif(NOT BOOST_INCLUDE_DIR)

# libraries, as specified in the COMPONENTS
foreach(component ${MyBoost_FIND_COMPONENTS})
  message( "Looking for Boost component ${component} as boost_${component}" )
  find_file( MyBoost_LIBRARY_${component}
    libboost_${component}.a
    PATH ${BOOST_DIR}/lib NO_DEFAULT_PATH)
  message(" MyBoost_LIBRARY_${component} : ${MyBoost_LIBRARY_${component}}")
  if(${MyBoost_LIBRARY_${component}} MATCHES "-NOTFOUND$")
    message(FATAL_ERROR "Didn't find boost ${component} library.")
  endif(${MyBoost_LIBRARY_${component}} MATCHES "-NOTFOUND$")
  set(BOOST_LIBRARIES ${BOOST_LIBRARIES} ${MyBoost_LIBRARY_${component}})
endforeach(component)

set(BOOST_FOUND "MyBoost_FOUND")

message("BOOST_INCLUDE_DIR : ${BOOST_INCLUDE_DIR}")
message("BOOST_LIBRARIES : ${BOOST_LIBRARIES}")
message("BOOST_FOUND : ${BOOST_FOUND}")
