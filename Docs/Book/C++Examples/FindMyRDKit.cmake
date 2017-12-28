# Attempts to find RDKit libraries using the current value of $RDBASE
# or, failing that, a version in my home directory
# It returns the static (.a) libraries not the .so ones because that's
# easiest for shipping (on Unix anyway. This may have to change once I start
# dealing with Windows as well.)
#
# It will define
# RDKIT_FOUND as MyRDKit_FOUND if it finds everything it needs
# RDKIT_INCLUDE_DIR
# RDKIT_LIBRARIES as requested

set(RDKIT_DIR $ENV{RDBASE})
if(NOT RDKIT_DIR)
  message( WARNING "Using RDKit at /home/cosgrove/RDKit_2013_09_1" )
  set(RDKIT_DIR "/home/cosgrove/RDKit_2013_09_1")
endif(NOT RDKIT_DIR)

set(RDKIT_INCLUDE_DIR ${RDKIT_DIR}/Code)

set(RDKIT_FOUND "MyRDKit_FOUND")
# libraries, as specified in the COMPONENTS
foreach(component ${MyRDKit_FIND_COMPONENTS})
  message( "Looking for RDKit component ${component}" )
  find_file( MyRDKit_LIBRARY_${component}
    libRDKit${component}.so
    PATH ${RDKIT_DIR}/lib NO_DEFAULT_PATH)
  message("MyRDKit_LIBRARY_${component} : ${MyRDKit_LIBRARY_${component}}")
  if(${MyRDKit_LIBRARY_${component}} MATCHES "-NOTFOUND$")
    message(FATAL_ERROR "Didn't find RDKit ${component} library.")
  endif(${MyRDKit_LIBRARY_${component}} MATCHES "-NOTFOUND$")
  set(RDKIT_LIBRARIES ${RDKIT_LIBRARIES} ${MyRDKit_LIBRARY_${component}})
endforeach(component)

message("RDKIT_INCLUDE_DIR : ${RDKIT_INCLUDE_DIR}")
message("RDKIT_LIBRARIES : ${RDKIT_LIBRARIES}")
message("RDKIT_FOUND : ${RDKIT_FOUND}")
