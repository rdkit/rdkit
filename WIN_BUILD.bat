call "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" amd64

rem CLEAN
rename build buildtodelete
rename lib libtodelete
del /f/s/q buildtodelete > nul
rmdir /s/q buildtodelete
del /f/s/q libtodelete > nul
rmdir /s/q libtodelete

rem CALL CMAKE
SET PYTHON_DIR=%NMRS_DEV_PYTHONBASE%
SET PYTHON_LIBRARY=%PYTHON_DIR%/libs/python35.lib
SET NUMPY_DIR=%PYTHON_DIR%/Lib/site-packages/numpy/core/include
SET BOOST_ROOT=%NMRS_DEV_BOOSTBASE%
SET BOOST_LIBRARYDIR=%BOOST_ROOT%/lib
mkdir build
cd build
cmake -DRDK_BUILD_PYTHON_WRAPPERS=ON -DPYTHON_EXECUTABLE=%PYTHON_DIR%/python.exe -DPYTHON_INCLUDE_DIR=%PYTHON_DIR%/include -DPYTHON_LIBRARY=%PYTHON_LIBRARY% -DPYTHON_NUMPY_INCLUDE_PATH=%NUMPY_DIR% -DBOOST_ROOT=%BOOST_ROOT% -DBOOST_LIBRARYDIR=%BOOST_LIBRARYDIR% -G"Visual Studio 14 2015 Win64" ..
cd ..

rem BUILD
cd build
MSBuild /m:12 /p:Configuration=Debug INSTALL.vcxproj
MSBuild /m:12 /p:Configuration=Release INSTALL.vcxproj
cd ..