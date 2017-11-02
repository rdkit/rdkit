set(RDKit_Year "2018")
set(RDKit_Month "03")
set(RDKit_Revision "1.dev1")
set(RDKit_ABI "1")

# Do not modify below this line
# if you need to, please check that pkg_version.py
# still extracts the correct data from this file
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac OS X specific code
  set(RDKit_VERSION "${RDKit_Year}.${RDKit_Month}")
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RDKit_VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(RDKit_RELEASENAME "${RDKit_Year}.${RDKit_Month}")
if (RDKit_Revision)
  set(RDKit_RELEASENAME "${RDKit_RELEASENAME}.${RDKit_Revision}")
  set(RDKit_VERSION "${RDKit_VERSION}.${RDKit_Revision}")
else(RDKit_Revision)
  set(RDKit_VERSION "${RDKit_VERSION}.0")
endif(RDKit_Revision)
