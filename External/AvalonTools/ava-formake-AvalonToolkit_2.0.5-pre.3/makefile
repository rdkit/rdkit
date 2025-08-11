SRCDIR := src/main
TESTDIR := src/test

JDK_INCLUDE ?= ${JAVA_HOME}/include

LD_OPTS := -fPIC -lm -fno-exceptions

# Detect OS flavors as required according to https://stackoverflow.com/questions/714100/os-detecting-makefile/14777895#14777895
ifeq ($(OS),Windows_NT)     # is Windows_NT on XP, 2000, 7, Vista, 10...
	# Windows flavor of operating systems
    FixPath = $(subst /,\,$1)	# replace slashes by backslashes. Might not be needed.
    LD_OPTS += -lole32
    DETECTED_OS := Windows
    # If path to javah is non-empty then it's safe to use it
	ifneq ($(realpath ${JAVA_HOME}/bin/javah.exe),) 
		JAVAH := TRUE
	endif
	KILL_AT = -Wl,--kill-at
	JDK_MACHINE_INCLUDE ?= ${JAVA_HOME}/include/win32
    MACHINE_DIR := win32
    GCC = i686-w64-mingw32-gcc
    ifeq ($(BITS),32)
        # Set cross compiler down to generate 32 bit artefacts
    	CC := $(GCC) -m32 -Wno-unused-but-set-variable -Wno-char-subscripts -Wno-misleading-indentation
    	LD := $(GCC) -m32 ${LD_OPTS} -static 
		WININCDIRS := ${SRCDIR}/C/include
    else
    	# Assuming the current architecture is 64 bit if BITS was anything other than 32
    	BITS := 64
    	# Patch to care for missing typedef of __int64
    	CC := $(GCC) -D__int64="long long" -Wno-unused-but-set-variable -Wno-char-subscripts -Wno-misleading-indentation
    	LD := $(GCC) -D__int64="long long" ${LD_OPTS} -static 
    endif
	LIB_PREFIX := lib
    EXE := .exe
    OBJ_EXT := obj
	SHARED_LIB_EXT := dll
	STATIC_LIB_EXT := a
	CLASS_PATH := "build/classes;lib/*"
	INCLUDEDIRS := ${SRCDIR}/C/include/
else 
    DETECTED_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
    # If path to javah is non-empty then it's safe to use it
	ifneq ($(realpath ${JAVA_HOME}/bin/javah),) 
		JAVAH := TRUE
	endif 
	ifeq ($(DETECTED_OS),Darwin)        # Mac OS X
		JDK_MACHINE_INCLUDE ?= ${JAVA_HOME}/include/darwin
	    MACHINE_DIR := darwin
	else
	    STATIC_LIB_FLAGS := -static-libgcc
	    DETECTED_OS := Linux
		JDK_MACHINE_INCLUDE ?= ${JAVA_HOME}/include/linux
	    MACHINE_DIR := linux
	endif
	LIB_PREFIX := lib
    ifeq ($(shell uname), Linux)
        RM = rm -f
        FixPath = $1
    endif
    ANT_HOME ?= /usr
    EXE := 
    OBJ_EXT := o
	SHARED_LIB_EXT := so
	STATIC_LIB_EXT := a
	CLASS_PATH := "build/classes:lib/*"
	INCLUDEDIRS := ${SRCDIR}/C/include
    CC := gcc -fPIC -fno-exceptions
    LD := gcc ${LD_OPTS} ${STATIC_LIB_FLAGS} 
endif

CSRCDIR := ${SRCDIR}/C/common
WINDIR := ${SRCDIR}/C/mswindows
LIBSRCDIR := ${SRCDIR}/C/library
PRGSRCDIR := ${SRCDIR}/C/programs
JNISRCDIR := ${SRCDIR}/C/jni
JAVASRCDIR := ${SRCDIR}/java

ANT := ${ANT_HOME}/bin/ant

include makefile.defs

BUILDDIR := build
LIBDIR := lib
TARGETDIR := target
EXETARGETS := ${TARGETDIR}/executables
LIBTARGETS := ${TARGETDIR}/libraries
JNIINCDIR := ${BUILDDIR}/jni_include

SMI2MOL_TARGETS := ${EXETARGETS}/smi2mol${EXE} ${EXETARGETS}/smi2rdf${EXE} ${EXETARGETS}/sma2mol${EXE} ${EXETARGETS}/mol2smi${EXE} ${EXETARGETS}/mol2sma${EXE} ${EXETARGETS}/mol2tbl${EXE}
#
PROGRAMS := ${EXETARGETS}/struchk${EXE} ${EXETARGETS}/canonizer${EXE} ${EXETARGETS}/matchtest${EXE} ${EXETARGETS}/test.exe ${SMI2MOL_TARGETS}

COMMON_OBJECTS := \
	${BUILDDIR}/forio.${OBJ_EXT} \
	${BUILDDIR}/geometry.${OBJ_EXT} \
	${BUILDDIR}/graph.${OBJ_EXT} \
	${BUILDDIR}/layout.${OBJ_EXT} \
	${BUILDDIR}/pattern.${OBJ_EXT} \
	${BUILDDIR}/perceive.${OBJ_EXT} \
	${BUILDDIR}/reaccsio.${OBJ_EXT} \
	${BUILDDIR}/set.${OBJ_EXT} \
	${BUILDDIR}/stereo.${OBJ_EXT} \
	${BUILDDIR}/symbol_lists.${OBJ_EXT} \
	${BUILDDIR}/symboltable.${OBJ_EXT} \
    ${BUILDDIR}/utilities.${OBJ_EXT} \
    ${BUILDDIR}/local.${OBJ_EXT}

RDKIT_OBJECTS := \
	${BUILDDIR}/layout.${OBJ_EXT} \
	${BUILDDIR}/symboltable.${OBJ_EXT} \
	${BUILDDIR}/patclean.${OBJ_EXT} \
    ${BUILDDIR}/utilities.${OBJ_EXT} \
	${BUILDDIR}/symbol_lists.${OBJ_EXT} \
	${BUILDDIR}/stereo.${OBJ_EXT} \
	${BUILDDIR}/set.${OBJ_EXT} \
	${BUILDDIR}/perceive.${OBJ_EXT} \
    ${BUILDDIR}/local.${OBJ_EXT} \
	${BUILDDIR}/graph.${OBJ_EXT} \
	${BUILDDIR}/geometry.${OBJ_EXT} \
	${BUILDDIR}/forio.${OBJ_EXT} \
	${BUILDDIR}/depictutil.${OBJ_EXT} \
	${BUILDDIR}/denormal.${OBJ_EXT} \
	${BUILDDIR}/casutils.${OBJ_EXT} \
	${BUILDDIR}/ssmatch.${OBJ_EXT} \
	${BUILDDIR}/rtutils.${OBJ_EXT} \
	${BUILDDIR}/smi2mol.${OBJ_EXT} \
	${BUILDDIR}/didepict.${OBJ_EXT} \
	${BUILDDIR}/pattern.${OBJ_EXT} \
	${BUILDDIR}/canonizer.${OBJ_EXT} \
	${BUILDDIR}/aacheck.${OBJ_EXT} \
	${BUILDDIR}/fixcharges.${OBJ_EXT} \
	${BUILDDIR}/struchk.${OBJ_EXT} \
	${BUILDDIR}/reaccsio.${OBJ_EXT} \
	${BUILDDIR}/hashcode.${OBJ_EXT}

# Build all targets that create artifacts. There are two further targets for Docker to
# start the depicter service directly. Note: The Docker system needs to be on the PATH
# for these additional targets to work.
all:	clean variables download programs javaclasses archives jni_libraries tomcat_lib
#
ifdef OS
readme :
	sed 's/\.exe/.exe/' <${SRCDIR}/resources/readme_platform.txt \
	   >${TARGETDIR}/readme_windows.txt
#
jni_libraries : ${LIBTARGETS}/toolkit${BITS}.${SHARED_LIB_EXT} ${LIBTARGETS}/depict${BITS}.${SHARED_LIB_EXT} \
                ${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} ${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} \
                ${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT} ${LIBTARGETS}/JNIWinTools.${SHARED_LIB_EXT} 
javaclasses :
	"${JAVA_HOME}/bin/javac" -cp ${CLASS_PATH} -d build/classes -sourcepath ${JAVASRCDIR} ${JAVASRCDIR}/avalon/jsp/servlets/*.java ${JAVASRCDIR}/jni/JNIWinTools.java
#
else
readme :
	sed 's/\.exe//' <${SRCDIR}/resources/readme_platform.txt \
	   >${TARGETDIR}/readme_linux.txt
#
jni_libraries : ${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} ${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} ${LIBTARGETS}/${LIB_PREFIX}rohde_clib.${SHARED_LIB_EXT} \
				${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT}
javaclasses :
	"${JAVA_HOME}/bin/javac" -cp ${CLASS_PATH} -d build/classes -sourcepath ${JAVASRCDIR} ${JAVASRCDIR}/avalon/jsp/servlets/*.java ${JAVASRCDIR}/jni/UnixTools.java
#
endif

#
variables: readme
	echo CLASS_PATH=${CLASS_PATH}
	echo OS=${OS}
	echo CYGWIN=${CYGWIN}
	echo BITS=${BITS}
	echo CC=${CC}
	echo LD=${LD}
	echo WINDIR=${WINDIR}
	echo DETECTED_OS=${DETECTED_OS}
	echo JAVAH=${JAVAH}
	echo 0: `which javac.exe`
	$(eval JAVAC=$(shell which javac.exe))
	echo 1: JAVAC=${JAVAC}

programs: ${PROGRAMS} ${LIBTARGETS}/avalon4rdkit.${STATIC_LIB_EXT} ${LIBTARGETS}/avalon_tools.${STATIC_LIB_EXT}

run_tests: programs
	${EXETARGETS}/canonizer${EXE} <${TESTDIR}/resources/canonizer.in >canonizer.out
	diff -Z ${TESTDIR}/resources/canonizer.out canonizer.out
#
ifdef JAVAH
javastubs : javastubs_javah
else
javastubs : javastubs_nojavah
endif
#
# Create header files for C code (pre-Java8)
javastubs_javah : download javaclasses 
	"${JAVA_HOME}/bin/javah" -classpath "${BUILDDIR}/classes" -d ${JNIINCDIR} jni.JNIDepict
	"${JAVA_HOME}/bin/javah" -classpath "${BUILDDIR}/classes" -d ${JNIINCDIR} jni.JNIMatch
	test ! -f "${BUILDDIR}/classes/jni/JNIWinTools.class" ||  "${JAVA_HOME}/bin/javah" -classpath "${BUILDDIR}/classes" -d ${JNIINCDIR} jni.JNIWinTools
	"${JAVA_HOME}/bin/javah" -classpath "${BUILDDIR}/classes" -d ${JNIINCDIR} avalon.jni.JNISmi2Mol
#
# Create header files for C code (post-Java8)
javastubs_nojavah : download javaclasses 
	"${JAVA_HOME}/bin/javac" -h ${JNIINCDIR} -d "${BUILDDIR}/classes" -cp "${BUILDDIR}/classes"  ${JAVASRCDIR}/jni/JNIDepict.java
	"${JAVA_HOME}/bin/javac" -h ${JNIINCDIR} -d "${BUILDDIR}/classes"  -cp "${BUILDDIR}/classes" ${JAVASRCDIR}/jni/JNIMatch.java
	test ! -f "${BUILDDIR}/classes/jni/JNIWinTools.class" ||  "${JAVA_HOME}/bin/javac" -h ${JNIINCDIR} -d "${BUILDDIR}/classes" -cp "${BUILDDIR}/classes"  ${JAVASRCDIR}/jni/JNIWinTools.java
	"${JAVA_HOME}/bin/javac" -h ${JNIINCDIR} -d "${BUILDDIR}/classes" -cp "${BUILDDIR}/classes" ${JAVASRCDIR}/avalon/jni/JNISmi2Mol.java
#
archives :
	cp ${SRCDIR}/resources/build.properties ${BUILDDIR}/classes
	${ANT_HOME}/bin/ant depictJar depictjniJar depictWar depictWarWithCORS
#
tomcat_lib : archives ${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} ${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} ${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT}
	cp ${TARGETDIR}/archives/depictjni.jar $(TARGETDIR)/tomcat_lib
	cp ${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} $(TARGETDIR)/tomcat_lib
	cp ${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} $(TARGETDIR)/tomcat_lib
	cp ${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT} $(TARGETDIR)/tomcat_lib
	# make sure Linux systems also find the file
	cp ${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT} $(TARGETDIR)/tomcat_lib/${LIB_PREFIX}JNISmi2Mol.${SHARED_LIB_EXT}
#
# Build the Docker image to be run by the tomcat_run target
tomcat_build : tomcat_lib
	docker build -f Dockerfile -t depicter:1.0.0 .
#
# Run the Tomcat Docker image serving depicter on the host's port 80
# The image is run detached, which means the Docker system needs to be used to stop it.
tomcat_run : tomcat_build
	docker run --rm -d --publish 80:8080 depicter:1.0.0
#
download :
	test -f ${LIBDIR}/log4j-1.2.8.jar || curl -o ${LIBDIR}/log4j-1.2.8.jar ${LOG4J_URL}
	test -f ${LIBDIR}/commons-logging-1.0.3.jar || curl -o ${LIBDIR}/commons-logging-1.0.3.jar ${COMMONS_URL}
	test -f ${LIBDIR}/cors-filter-1.3.2.jar || curl -o ${LIBDIR}/cors-filter-1.3.2.jar ${CORS_URL}
	test -f ${LIBDIR}/java-property-utils-1.6.jar || curl -o ${LIBDIR}/java-property-utils-1.6.jar ${PROP_UTILS_URL}
	test -f ${LIBDIR}/servlet-api-2.4.jar || curl -o ${LIBDIR}/servlet-api-2.4.jar ${SERVLET_URL}
	test -f ${LIBDIR}/jstl-1.2.jar || curl -o ${LIBDIR}/jstl-1.2.jar ${JSTL_URL}
#
clean:
	@echo " Cleaning..."
	@echo " $(RM) -r $(BUILDDIR) $(TARGETDIR) $(LIBDIR)"; $(RM) -r $(BUILDDIR) $(TARGETDIR) $(LIBDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(BUILDDIR)/classes
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(TARGETDIR)/tomcat_lib
	@mkdir -p $(EXETARGETS)
	@mkdir -p $(LIBTARGETS)
#
# Individual object files
${BUILDDIR}/struchk.${OBJ_EXT}:	${PRGSRCDIR}/struchk.c
	${CC} -I${INCLUDEDIRS} -c -o $@ $<
#
# Catch-all object files
${BUILDDIR}/%.${OBJ_EXT} : ${CSRCDIR}/%.c
	${CC} -I${INCLUDEDIRS} -c $(CFLAGS) $< -o $@
#
${EXETARGETS}/struchk${EXE} : ${PRGSRCDIR}/struchk.c ${COMMON_OBJECTS} \
	${BUILDDIR}/hashcode.${OBJ_EXT} \
	${BUILDDIR}/patclean.${OBJ_EXT} \
	${BUILDDIR}/aacheck.${OBJ_EXT} \
	${BUILDDIR}/fixcharges.${OBJ_EXT} \
    ${BUILDDIR}/ssmatch.${OBJ_EXT}
	${CC} -DMAIN -I${INCLUDEDIRS} -o ${EXETARGETS}/struchk${EXE} \
	   ${PRGSRCDIR}/struchk.c ${COMMON_OBJECTS} \
	   ${BUILDDIR}/hashcode.${OBJ_EXT} \
	   ${BUILDDIR}/patclean.${OBJ_EXT} \
	   ${BUILDDIR}/aacheck.${OBJ_EXT} \
	   ${BUILDDIR}/fixcharges.${OBJ_EXT} \
	   ${BUILDDIR}/ssmatch.${OBJ_EXT} ${LD_OPTS}
#
${LIBTARGETS}/avalon4rdkit.${STATIC_LIB_EXT} : ${RDKIT_OBJECTS}
	ar rcs ${LIBTARGETS}/libavalon4rdkit.${STATIC_LIB_EXT} ${RDKIT_OBJECTS}
#
${LIBTARGETS}/avalon_tools.${STATIC_LIB_EXT} : ${BUILDDIR}/*.${OBJ_EXT}
	ar rcs ${LIBTARGETS}/libavalon_tools.${STATIC_LIB_EXT} ${BUILDDIR}/*.${OBJ_EXT}
#
${EXETARGETS}/matchtest${EXE} : ${PRGSRCDIR}/matchtest.c ${COMMON_OBJECTS} \
	${BUILDDIR}/hashcode.${OBJ_EXT} \
	${BUILDDIR}/patclean.${OBJ_EXT} \
	${BUILDDIR}/fixcharges.${OBJ_EXT} \
    ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	${BUILDDIR}/casutils.${OBJ_EXT} \
	${BUILDDIR}/denormal.${OBJ_EXT} \
	${BUILDDIR}/rtutils.${OBJ_EXT} \
	${BUILDDIR}/smi2mol.${OBJ_EXT} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/matchtest${EXE} \
	   ${PRGSRCDIR}/matchtest.c ${COMMON_OBJECTS} \
	   ${BUILDDIR}/hashcode.${OBJ_EXT} \
	   ${BUILDDIR}/patclean.${OBJ_EXT} \
	   ${BUILDDIR}/fixcharges.${OBJ_EXT} \
	   ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	   ${BUILDDIR}/casutils.${OBJ_EXT} \
	   ${BUILDDIR}/denormal.${OBJ_EXT} \
	   ${BUILDDIR}/rtutils.${OBJ_EXT} \
	   ${BUILDDIR}/smi2mol.${OBJ_EXT} ${LD_OPTS}
#
${EXETARGETS}/canonizer${EXE} : ${PRGSRCDIR}/canonizer_main.c ${COMMON_OBJECTS} \
	${BUILDDIR}/canonizer.${OBJ_EXT} \
	${BUILDDIR}/casutils.${OBJ_EXT} \
	${BUILDDIR}/denormal.${OBJ_EXT} \
	${BUILDDIR}/rtutils.${OBJ_EXT} \
	${BUILDDIR}/smi2mol.${OBJ_EXT} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/canonizer${EXE} \
	   ${PRGSRCDIR}/canonizer_main.c ${COMMON_OBJECTS} \
	   ${BUILDDIR}/canonizer.${OBJ_EXT} \
	   ${BUILDDIR}/casutils.${OBJ_EXT} \
	   ${BUILDDIR}/denormal.${OBJ_EXT} \
	   ${BUILDDIR}/rtutils.${OBJ_EXT} \
	   ${BUILDDIR}/smi2mol.${OBJ_EXT} ${LD_OPTS}
#
SMI2MOL_OBJECTS := \
	${BUILDDIR}/canonizer.${OBJ_EXT} \
	${BUILDDIR}/casutils.${OBJ_EXT} \
	${BUILDDIR}/denormal.${OBJ_EXT} \
	${BUILDDIR}/ssmatch.${OBJ_EXT} \
	${BUILDDIR}/didepict.${OBJ_EXT} \
	${BUILDDIR}/shortcut.${OBJ_EXT} \
	${BUILDDIR}/hashcode.${OBJ_EXT} \
	${BUILDDIR}/rtutils.${OBJ_EXT} \
	${BUILDDIR}/smi2mol.${OBJ_EXT} 
#
${EXETARGETS}/smi2mol${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/smi2mol${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/sma2mol${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/sma2mol${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/mol2sma${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/mol2sma${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/smi2plt${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/smi2plt${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/smi2rdf${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/smi2rdf${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/smi2smi${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/smi2smi${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/mol2tbl${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/mol2tbl${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/rdf2smi${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/rdf2smi${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/mol2smi${EXE} : ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} 
	${CC} -I${INCLUDEDIRS} -o ${EXETARGETS}/mol2smi${EXE} \
	   ${PRGSRCDIR}/smi2mol_main.c ${COMMON_OBJECTS} ${SMI2MOL_OBJECTS} ${LD_OPTS}
#
${EXETARGETS}/test.exe: ${BUILDDIR}/test.obj
	${CC} -o ${EXETARGETS}/test.exe ${BUILDDIR}/test.obj
#
${BUILDDIR}/test.obj:	${TESTDIR}/C/test.c
	${CC} -c -o $@ $<

# Individual object files
${BUILDDIR}/toolkit${BITS}.${OBJ_EXT} : ${WINDIR}/toolkit.c  ${WINDIR}/toolkit.h
	${CC} -c -I ${WINDIR} \
	       --verbose -I ${INCLUDEDIRS} \
       	       -o ${BUILDDIR}/toolkit${BITS}.${OBJ_EXT} ${WINDIR}/toolkit.c
#
${BUILDDIR}/depict${BITS}.${OBJ_EXT} : ${WINDIR}/depict.c  ${WINDIR}/depict.h
	${CC} -c -I ${WINDIR} \
	       --verbose -I ${INCLUDEDIRS} \
			-o ${BUILDDIR}/depict${BITS}.${OBJ_EXT} ${WINDIR}/depict.c
#
# This is a useful URL describing among other calling convention and name mangling issues the kill-at option
# http://www.willus.com/mingw/yongweiwu_stdcall.html
${LIBTARGETS}/depict${BITS}.${SHARED_LIB_EXT} : \
    	${COMMON_OBJECTS} \
    	${SMI2MOL_OBJECTS} \
    	${BUILDDIR}/sketch.${OBJ_EXT} \
    	${BUILDDIR}/depict${BITS}.${OBJ_EXT}
	${LD} -shared -o ${LIBTARGETS}/depict${BITS}.${SHARED_LIB_EXT} \
	--verbose \
    ${KILL_AT} \
    -mwindows \
    ${COMMON_OBJECTS} \
    ${SMI2MOL_OBJECTS} ${BUILDDIR}/sketch.${OBJ_EXT} \
    ${BUILDDIR}/depict${BITS}.${OBJ_EXT} \
	${LD_OPTS}
#
${LIBTARGETS}/toolkit${BITS}.${SHARED_LIB_EXT} : \
    	${COMMON_OBJECTS} \
    	${SMI2MOL_OBJECTS} \
    	${BUILDDIR}/toolkit${BITS}.${OBJ_EXT}
	${LD} -shared -o ${LIBTARGETS}/toolkit${BITS}.${SHARED_LIB_EXT} \
	--verbose \
    ${KILL_AT} \
    -mwindows \
    ${COMMON_OBJECTS} \
    ${SMI2MOL_OBJECTS} \
    ${BUILDDIR}/toolkit${BITS}.${OBJ_EXT} \
	${LD_OPTS}
#
# JNI section
${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT} : \
	       ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/depictutil.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} ${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/hashcode.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} ${BUILDDIR}/patclean.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} ${BUILDDIR}/reaccsio.${OBJ_EXT} \
	       ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
	       ${BUILDDIR}/shortcut.${OBJ_EXT} \
	       ${JNISRCDIR}/avalon_jni_JNISmi2Mol.cpp  javastubs
	${LD} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I "${JDK_INCLUDE}" \
           -I "${JDK_MACHINE_INCLUDE}" \
           -I "${JNIINCDIR}" \
	       -I ${INCLUDEDIRS} \
	      -shared ${JNISRCDIR}/avalon_jni_JNISmi2Mol.cpp \
	       ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/depictutil.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} ${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/hashcode.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} ${BUILDDIR}/patclean.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} ${BUILDDIR}/reaccsio.${OBJ_EXT} \
	       ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
	       ${BUILDDIR}/shortcut.${OBJ_EXT} \
	      -o ${LIBTARGETS}/avalon_jni_JNISmi2Mol.${SHARED_LIB_EXT} \
		   ${LD_OPTS}
#
${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} : \
	       ${BUILDDIR}/canonizer.${OBJ_EXT} ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} ${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/hashcode.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} ${BUILDDIR}/reaccsio.${OBJ_EXT} \
	       ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
	      ${JNISRCDIR}/jni_JNIMatch.cpp  javastubs
	${LD} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I "${JDK_INCLUDE}" \
           -I "${JDK_MACHINE_INCLUDE}" \
           -I "${JNIINCDIR}" \
	       -I ${INCLUDEDIRS} \
              -shared ${JNISRCDIR}/jni_JNIMatch.cpp \
	       ${BUILDDIR}/canonizer.${OBJ_EXT} ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} ${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/hashcode.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} ${BUILDDIR}/reaccsio.${OBJ_EXT} \
	       ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
		   -o ${LIBTARGETS}/${LIB_PREFIX}JNIMatch.${SHARED_LIB_EXT} \
		   ${LD_OPTS}
#
${LIBTARGETS}/JNIWinTools.${SHARED_LIB_EXT} : \
	       ${BUILDDIR}/depictutil.${OBJ_EXT} \
	      ${JNISRCDIR}/jni_JNIWinTools.cpp  javastubs
	${CC} -Wall -mwindows -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I "${JDK_INCLUDE}" \
           -I "${JDK_MACHINE_INCLUDE}" \
           -I "${JNIINCDIR}" \
		   -I ${WINDIR} \
	       -I ${INCLUDEDIRS} \
	       ${BUILDDIR}/casutils.${OBJ_EXT} \
	       ${BUILDDIR}/depictutil.${OBJ_EXT} \
	       ${BUILDDIR}/denormal.${OBJ_EXT} \
	       ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} \
	       ${BUILDDIR}/graph.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} \
	       ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/perceive.${OBJ_EXT} \
	       ${BUILDDIR}/reaccsio.${OBJ_EXT} \
	       ${BUILDDIR}/rtutils.${OBJ_EXT} \
	       ${BUILDDIR}/set.${OBJ_EXT} \
	       ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} \
	       ${BUILDDIR}/symbol_lists.${OBJ_EXT} \
	       ${BUILDDIR}/utilities.${OBJ_EXT} \
		   -shared ${JNISRCDIR}/jni_JNIWinTools.cpp \
		   -o ${LIBTARGETS}/JNIWinTools.${SHARED_LIB_EXT} \
		   ${LD_OPTS}
#
${LIBTARGETS}/ClipboardTools.${SHARED_LIB_EXT} : \
	      ${JNISRCDIR}/jni_ClipboardTools.cpp  javastubs
	${CC} -Wall -mwindows -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I "${JDK_INCLUDE}" \
          -I "${JDK_MACHINE_INCLUDE}" \
	       -I ${INCLUDEDIRS} \
           -shared ${JNISRCDIR}/jni_ClipboardTools.cpp \
		   -o ${LIBTARGETS}/ClipboardTools.${SHARED_LIB_EXT} \
		   ${LD_OPTS}
#
${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} : \
	       ${BUILDDIR}/casutils.${OBJ_EXT} \
	       ${BUILDDIR}/denormal.${OBJ_EXT} \
		   ${BUILDDIR}/depictutil.${OBJ_EXT} \
		   ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} \
		   ${BUILDDIR}/graph.${OBJ_EXT} \
		   ${BUILDDIR}/hashcode.${OBJ_EXT} \
		   ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} \
		   ${BUILDDIR}/patclean.${OBJ_EXT} \
		   ${BUILDDIR}/perceive.${OBJ_EXT} \
		   ${BUILDDIR}/reaccsio.${OBJ_EXT} \
		   ${BUILDDIR}/rtutils.${OBJ_EXT} \
	       ${BUILDDIR}/set.${OBJ_EXT} \
		   ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} \
		   ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} \
		   ${BUILDDIR}/symbol_lists.${OBJ_EXT} \
		   ${BUILDDIR}/symboltable.${OBJ_EXT} \
		   ${BUILDDIR}/utilities.${OBJ_EXT} \
	      ${JNISRCDIR}/jni_JNIDepict.cpp  javastubs
	${CC} -Wall -D_JNI_IMPLEMENTATION_ ${KILL_AT} \
	      -I "${JDK_INCLUDE}" \
          -I "${JDK_MACHINE_INCLUDE}" \
	       -I ${JNIINCDIR} \
	       -I ${INCLUDEDIRS} \
              -shared ${JNISRCDIR}/jni_JNIDepict.cpp \
	       ${BUILDDIR}/casutils.${OBJ_EXT} \
	       ${BUILDDIR}/denormal.${OBJ_EXT} \
		   ${BUILDDIR}/depictutil.${OBJ_EXT} \
		   ${BUILDDIR}/forio.${OBJ_EXT} \
	       ${BUILDDIR}/geometry.${OBJ_EXT} \
		   ${BUILDDIR}/graph.${OBJ_EXT} \
		   ${BUILDDIR}/hashcode.${OBJ_EXT} \
		   ${BUILDDIR}/layout.${OBJ_EXT} \
	       ${BUILDDIR}/local.${OBJ_EXT} \
		   ${BUILDDIR}/patclean.${OBJ_EXT} \
		   ${BUILDDIR}/perceive.${OBJ_EXT} \
		   ${BUILDDIR}/reaccsio.${OBJ_EXT} \
		   ${BUILDDIR}/rtutils.${OBJ_EXT} \
	       ${BUILDDIR}/set.${OBJ_EXT} \
		   ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} \
		   ${BUILDDIR}/ssmatch.${OBJ_EXT} \
	       ${BUILDDIR}/stereo.${OBJ_EXT} \
		   ${BUILDDIR}/symbol_lists.${OBJ_EXT} \
		   ${BUILDDIR}/symboltable.${OBJ_EXT} \
		   ${BUILDDIR}/utilities.${OBJ_EXT} \
		   -o ${LIBTARGETS}/${LIB_PREFIX}JNIDepict.${SHARED_LIB_EXT} \
		   ${LD_OPTS}
#
${LIBTARGETS}/librohde_clib.${SHARED_LIB_EXT} : \
	${BUILDDIR}/aacheck.${OBJ_EXT} \
	${BUILDDIR}/fixcharges.${OBJ_EXT} \
	${BUILDDIR}/hashcode.${OBJ_EXT} \
	${BUILDDIR}/patclean.${OBJ_EXT} \
	${BUILDDIR}/ssmatch.${OBJ_EXT} \
	${BUILDDIR}/canonizer.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} ${BUILDDIR}/geometry.${OBJ_EXT} \
	${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} \
	${BUILDDIR}/reaccsio.${OBJ_EXT} ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} \
	${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
	${BUILDDIR}/local.${OBJ_EXT}
	${CC} -Wall \
	       -I ${INCLUDEDIRS} \
	      -o ${LIBTARGETS}/librohde_clib.${SHARED_LIB_EXT} \
	      -shared ${PRGSRCDIR}/struchk.c \
		${BUILDDIR}/aacheck.${OBJ_EXT} ${BUILDDIR}/canonizer.${OBJ_EXT} ${BUILDDIR}/smi2mol.${OBJ_EXT} ${BUILDDIR}/pattern.${OBJ_EXT} ${BUILDDIR}/denormal.${OBJ_EXT} ${BUILDDIR}/casutils.${OBJ_EXT} ${BUILDDIR}/fixcharges.${OBJ_EXT} ${BUILDDIR}/forio.${OBJ_EXT} ${BUILDDIR}/geometry.${OBJ_EXT} \
		${BUILDDIR}/graph.${OBJ_EXT} ${BUILDDIR}/hashcode.${OBJ_EXT} ${BUILDDIR}/layout.${OBJ_EXT} \
		${BUILDDIR}/patclean.${OBJ_EXT} ${BUILDDIR}/perceive.${OBJ_EXT} ${BUILDDIR}/reaccsio.${OBJ_EXT} ${BUILDDIR}/rtutils.${OBJ_EXT} ${BUILDDIR}/set.${OBJ_EXT} \
		${BUILDDIR}/ssmatch.${OBJ_EXT} ${BUILDDIR}/stereo.${OBJ_EXT} ${BUILDDIR}/symbol_lists.${OBJ_EXT} ${BUILDDIR}/symboltable.${OBJ_EXT} ${BUILDDIR}/utilities.${OBJ_EXT} \
	    ${BUILDDIR}/local.${OBJ_EXT} \
	    ${LD_OPTS}
