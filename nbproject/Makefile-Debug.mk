#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/CFastqCreator.o \
	${OBJECTDIR}/CIsbtGt.o \
	${OBJECTDIR}/CIsbtGt2Pt.o \
	${OBJECTDIR}/CIsbtGt2PtHit.o \
	${OBJECTDIR}/CIsbtGtAllele.o \
	${OBJECTDIR}/CIsbtPtAllele.o \
	${OBJECTDIR}/CIsbtVariant.o \
	${OBJECTDIR}/CMakeTrainingVcf.o \
	${OBJECTDIR}/CVariantChain.o \
	${OBJECTDIR}/CVariantChainVariation.o \
	${OBJECTDIR}/CVariantChains.o \
	${OBJECTDIR}/CVcf.o \
	${OBJECTDIR}/CVcfSnp.o \
	${OBJECTDIR}/ISBTAnno.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../bamtools/build/src/api/ -L../htslib -Wl,-rpath,'../bamtools/build/src/api/' -Wl,-rpath,'../htslib' ../bamtools/build/src/api/libbamtools.a -Wl,-rpath,'../MyTools/dist/Debug/GNU-Linux' -L../MyTools/dist/Debug/GNU-Linux -lMyTools

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood: ../bamtools/build/src/api/libbamtools.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood: ../MyTools/dist/Debug/GNU-Linux/libMyTools.so

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood ${OBJECTFILES} ${LDLIBSOPTIONS} -lz -lhts

${OBJECTDIR}/CFastqCreator.o: CFastqCreator.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CFastqCreator.o CFastqCreator.cpp

${OBJECTDIR}/CIsbtGt.o: CIsbtGt.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt.o CIsbtGt.cpp

${OBJECTDIR}/CIsbtGt2Pt.o: CIsbtGt2Pt.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt2Pt.o CIsbtGt2Pt.cpp

${OBJECTDIR}/CIsbtGt2PtHit.o: CIsbtGt2PtHit.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt2PtHit.o CIsbtGt2PtHit.cpp

${OBJECTDIR}/CIsbtGtAllele.o: CIsbtGtAllele.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGtAllele.o CIsbtGtAllele.cpp

${OBJECTDIR}/CIsbtPtAllele.o: CIsbtPtAllele.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtPtAllele.o CIsbtPtAllele.cpp

${OBJECTDIR}/CIsbtVariant.o: CIsbtVariant.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtVariant.o CIsbtVariant.cpp

${OBJECTDIR}/CMakeTrainingVcf.o: CMakeTrainingVcf.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CMakeTrainingVcf.o CMakeTrainingVcf.cpp

${OBJECTDIR}/CVariantChain.o: CVariantChain.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChain.o CVariantChain.cpp

${OBJECTDIR}/CVariantChainVariation.o: CVariantChainVariation.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChainVariation.o CVariantChainVariation.cpp

${OBJECTDIR}/CVariantChains.o: CVariantChains.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChains.o CVariantChains.cpp

${OBJECTDIR}/CVcf.o: CVcf.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVcf.o CVcf.cpp

${OBJECTDIR}/CVcfSnp.o: CVcfSnp.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVcfSnp.o CVcfSnp.cpp

${OBJECTDIR}/ISBTAnno.o: ISBTAnno.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ISBTAnno.o ISBTAnno.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../MyTools -I../htslib/htslib -I../htslib -I../bamtools/src -I../tclap-1.2.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:
	cd ../MyTools && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libMyTools.so
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood

# Subprojects
.clean-subprojects:
	cd ../MyTools && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
