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
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/CBigWigReader.o \
	${OBJECTDIR}/CFastqCreator.o \
	${OBJECTDIR}/CFastqReader.o \
	${OBJECTDIR}/CIsbtGt.o \
	${OBJECTDIR}/CIsbtGt2Pt.o \
	${OBJECTDIR}/CIsbtGt2PtHit.o \
	${OBJECTDIR}/CIsbtGtAllele.o \
	${OBJECTDIR}/CIsbtPtAllele.o \
	${OBJECTDIR}/CIsbtVariant.o \
	${OBJECTDIR}/CMakeTrainingVcf.o \
	${OBJECTDIR}/CMotifFinder.o \
	${OBJECTDIR}/CScoreHaplotype.o \
	${OBJECTDIR}/CTranscript.o \
	${OBJECTDIR}/CTranscriptAnno.o \
	${OBJECTDIR}/CVariantChain.o \
	${OBJECTDIR}/CVariantChainVariation.o \
	${OBJECTDIR}/CVariantChains.o \
	${OBJECTDIR}/CVcf.o \
	${OBJECTDIR}/CVcfSnp.o \
	${OBJECTDIR}/ISBTAnno.o \
	${OBJECTDIR}/gzstream.o \
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
LDLIBSOPTIONS=-L../htslib -L../libBigWig -Wl,-rpath,'../htslib' -Wl,-rpath,'../libBigWig' ../MyTools/dist/Release/GNU-Linux/libmytools.a -lhts -lBigWig

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood: ../MyTools/dist/Release/GNU-Linux/libmytools.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/deepblood ${OBJECTFILES} ${LDLIBSOPTIONS} -lz -pthread

${OBJECTDIR}/CBigWigReader.o: CBigWigReader.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CBigWigReader.o CBigWigReader.cpp

${OBJECTDIR}/CFastqCreator.o: CFastqCreator.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CFastqCreator.o CFastqCreator.cpp

${OBJECTDIR}/CFastqReader.o: CFastqReader.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CFastqReader.o CFastqReader.cpp

${OBJECTDIR}/CIsbtGt.o: CIsbtGt.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt.o CIsbtGt.cpp

${OBJECTDIR}/CIsbtGt2Pt.o: CIsbtGt2Pt.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt2Pt.o CIsbtGt2Pt.cpp

${OBJECTDIR}/CIsbtGt2PtHit.o: CIsbtGt2PtHit.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGt2PtHit.o CIsbtGt2PtHit.cpp

${OBJECTDIR}/CIsbtGtAllele.o: CIsbtGtAllele.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtGtAllele.o CIsbtGtAllele.cpp

${OBJECTDIR}/CIsbtPtAllele.o: CIsbtPtAllele.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtPtAllele.o CIsbtPtAllele.cpp

${OBJECTDIR}/CIsbtVariant.o: CIsbtVariant.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CIsbtVariant.o CIsbtVariant.cpp

${OBJECTDIR}/CMakeTrainingVcf.o: CMakeTrainingVcf.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CMakeTrainingVcf.o CMakeTrainingVcf.cpp

${OBJECTDIR}/CMotifFinder.o: CMotifFinder.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CMotifFinder.o CMotifFinder.cpp

${OBJECTDIR}/CScoreHaplotype.o: CScoreHaplotype.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CScoreHaplotype.o CScoreHaplotype.cpp

${OBJECTDIR}/CTranscript.o: CTranscript.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CTranscript.o CTranscript.cpp

${OBJECTDIR}/CTranscriptAnno.o: CTranscriptAnno.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CTranscriptAnno.o CTranscriptAnno.cpp

${OBJECTDIR}/CVariantChain.o: CVariantChain.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChain.o CVariantChain.cpp

${OBJECTDIR}/CVariantChainVariation.o: CVariantChainVariation.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChainVariation.o CVariantChainVariation.cpp

${OBJECTDIR}/CVariantChains.o: CVariantChains.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVariantChains.o CVariantChains.cpp

${OBJECTDIR}/CVcf.o: CVcf.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVcf.o CVcf.cpp

${OBJECTDIR}/CVcfSnp.o: CVcfSnp.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CVcfSnp.o CVcfSnp.cpp

${OBJECTDIR}/ISBTAnno.o: ISBTAnno.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ISBTAnno.o ISBTAnno.cpp

${OBJECTDIR}/gzstream.o: gzstream.C
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/gzstream.o gzstream.C

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I../MyTools -I../htslib/htslib -I../htslib -I../libBigWig -Itclap/include -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:
	cd ../MyTools && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:
	cd ../MyTools && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
