

# 
# Uncomment the following line if your compiler lacks the bool type
#
# -DNEED_BOOL


# 
# Optimization flags
#
OPTFLAGS= -O6
#OPTFLAGS=

#
# Debugging flags. -g adds debugging info
# Defining NDEBUG turns assertions off.
# -pg adds profiling info 
# 
#DEBUGFLAGS= -g 
# DEBUGFLAGS=
# DEBUGFLAGS= -DNDEBUG=1
# DEBUGFLAGS= -g -pg


#
# Under g++, -Wall enables all the warning messages
#
OTHERFLAGS= -Wall

#---------------------------------------------------------------
# Don't edit below this line
#---------------------------------------------------------------

LIBRARY=libvf.a

CXXFLAGS= -Iinclude $(OPTFLAGS) $(DEBUGFLAGS) $(OTHERFLAGS)

OBJS=	src/argedit.o src/argraph.o src/argloader.o src/error.o src/gene.o \
	src/gene_mesh.o src/match.o src/ull_state.o src/ull_sub_state.o \
	src/vf_mono_state.o src/vf_state.o src/vf_sub_state.o \
	src/vf2_state.o src/vf2_sub_state.o src/vf2_mono_state.o \
	src/sd_state.o \
	src/sortnodes.o src/xsubgraph.o

all:	lib/$(LIBRARY)

lib/$(LIBRARY):	$(OBJS)
	-rm lib/$(LIBRARY)
	ar rv lib/$(LIBRARY) $(OBJS)
	ranlib lib/$(LIBRARY)


tgz:
	-rm vflib2.tgz
	tar cvfz vflib2.tgz README Makefile src/*.cc src/*.h include/* doc/* lib/dummy

depend:
	makedepend -Iinclude -Y src/*

clean:
	-rm src/*.o

# DO NOT DELETE

src/argedit.o: include/argedit.h include/argraph.h include/error.h
src/argloader.o: include/argloader.h include/argedit.h include/argraph.h
src/argloader.o: include/allocpool.h include/error.h
src/argraph.o: include/argraph.h include/error.h
src/error.o: include/error.h
src/gene.o: include/argraph.h include/argedit.h include/error.h
src/gene.o: include/gene.h
src/gene_mesh.o: include/argraph.h include/argedit.h include/error.h
src/gene_mesh.o: include/gene.h
src/match.o: include/argraph.h include/match.h include/state.h
src/match.o: include/error.h
src/sd_state.o: include/argraph.h include/sd_state.h include/state.h
src/sd_state.cc.o: include/argraph.h include/sd_state.h include/state.h
src/sortnodes.o: src/sortnodes.h include/argraph.h
src/sortnodes.o: include/argraph.h
src/ull_state.o: include/ull_state.h include/argraph.h include/state.h
src/ull_state.o: include/error.h
src/ull_sub_state.o: include/ull_sub_state.h include/argraph.h
src/ull_sub_state.o: include/state.h include/error.h
src/vf2_mono_state.o: include/vf2_mono_state.h include/argraph.h
src/vf2_mono_state.o: include/state.h src/sortnodes.h include/error.h
src/vf2_state.o: include/vf2_state.h include/argraph.h include/state.h
src/vf2_state.o: include/error.h src/sortnodes.h
src/vf2_state.cc.o: include/vf2_state.h include/argraph.h include/state.h
src/vf2_state.cc.o: include/error.h src/sortnodes.h
src/vf2_sub_state.o: include/vf2_sub_state.h include/argraph.h
src/vf2_sub_state.o: include/state.h src/sortnodes.h include/error.h
src/vf_mono_state.o: include/vf_mono_state.h include/argraph.h
src/vf_mono_state.o: include/state.h include/error.h
src/vf_state.o: include/vf_state.h include/argraph.h include/state.h
src/vf_state.o: include/error.h
src/vf_sub_state.o: include/vf_sub_state.h include/argraph.h include/state.h
src/vf_sub_state.o: include/error.h
src/xsubgraph.o: include/argraph.h include/argedit.h include/xsubgraph.h
src/xsubgraph.o: include/error.h
