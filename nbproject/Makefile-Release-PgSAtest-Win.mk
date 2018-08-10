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
CND_PLATFORM=MinGW_1-Windows
CND_DLIB_EXT=dll
CND_CONF=Release-PgSAtest-Win
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include MakefileNetBeans

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/PgSAtest.o \
	${OBJECTDIR}/src/helper.o \
	${OBJECTDIR}/src/index/PgSAIndexFactory.o \
	${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o \
	${OBJECTDIR}/src/pgsaconfig.o \
	${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o \
	${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o \
	${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o \
	${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o \
	${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o \
	${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o \
	${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o \
	${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o \
	${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o \
	${OBJECTDIR}/src/readsset/DefaultReadsSet.o \
	${OBJECTDIR}/src/readsset/PackedReadsSet.o \
	${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o \
	${OBJECTDIR}/src/sais/sais.o \
	${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o \
	${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o \
	${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o \
	${OBJECTDIR}/src/test/PgSAtests.o \
	${OBJECTDIR}/src/test/testdata.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f3

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/src/tests/PgSASearchTest.o \
	${TESTDIR}/src/tests/PgSAgen-fasta-human.o

# C Compiler Flags
CFLAGS=-fomit-frame-pointer -Wall

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/PgSAtest.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/PgSAtest.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/PgSAtest ${OBJECTFILES} ${LDLIBSOPTIONS} -s

${OBJECTDIR}/src/PgSAtest.o: src/PgSAtest.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PgSAtest.o src/PgSAtest.cpp

${OBJECTDIR}/src/helper.o: src/helper.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/helper.o src/helper.cpp

${OBJECTDIR}/src/index/PgSAIndexFactory.o: src/index/PgSAIndexFactory.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/index
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/index/PgSAIndexFactory.o src/index/PgSAIndexFactory.cpp

${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o: src/index/cache/persistence/CountQueriesCachePersistence.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/index/cache/persistence
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o src/index/cache/persistence/CountQueriesCachePersistence.cpp

${OBJECTDIR}/src/pgsaconfig.o: src/pgsaconfig.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pgsaconfig.o src/pgsaconfig.cpp

${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o: src/pseudogenome/DefaultPseudoGenome.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o src/pseudogenome/DefaultPseudoGenome.cpp

${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o: src/pseudogenome/PackedPseudoGenome.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o src/pseudogenome/PackedPseudoGenome.cpp

${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o: src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp

${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o: src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.cpp

${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o: src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp

${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o: src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp

${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o: src/pseudogenome/packing/SymbolsPackingFacility.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/packing
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o src/pseudogenome/packing/SymbolsPackingFacility.cpp

${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o: src/pseudogenome/persistence/PseudoGenomePersistence.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/persistence
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o src/pseudogenome/persistence/PseudoGenomePersistence.cpp

${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o: src/pseudogenome/readslist/ListOfConstantLengthReads.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/readslist
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o src/pseudogenome/readslist/ListOfConstantLengthReads.cpp

${OBJECTDIR}/src/readsset/DefaultReadsSet.o: src/readsset/DefaultReadsSet.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/readsset
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/DefaultReadsSet.o src/readsset/DefaultReadsSet.cpp

${OBJECTDIR}/src/readsset/PackedReadsSet.o: src/readsset/PackedReadsSet.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/readsset
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/PackedReadsSet.o src/readsset/PackedReadsSet.cpp

${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o: src/readsset/iterator/ReadsSetIterator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/readsset/iterator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o src/readsset/iterator/ReadsSetIterator.cpp

${OBJECTDIR}/src/sais/sais.o: src/sais/sais.c nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/sais
	${RM} "$@.d"
	$(COMPILE.c) -I. -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sais/sais.o src/sais/sais.c

${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o: src/suffixarray/DefaultSuffixArray.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o src/suffixarray/DefaultSuffixArray.cpp

${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o: src/suffixarray/SparseSuffixArray.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o src/suffixarray/SparseSuffixArray.cpp

${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o: src/suffixarray/persistence/SuffixArrayPersistence.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray/persistence
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o src/suffixarray/persistence/SuffixArrayPersistence.cpp

${OBJECTDIR}/src/test/PgSAtests.o: src/test/PgSAtests.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/test
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/test/PgSAtests.o src/test/PgSAtests.cpp

${OBJECTDIR}/src/test/testdata.o: src/test/testdata.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/src/test
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/test/testdata.o src/test/testdata.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/f1: ${TESTDIR}/src/tests/PgSAgen-fasta-human.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS}   

${TESTDIR}/TestFiles/f3: ${TESTDIR}/src/tests/PgSASearchTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS}   


${TESTDIR}/src/tests/PgSAgen-fasta-human.o: src/tests/PgSAgen-fasta-human.cpp 
	${MKDIR} -p ${TESTDIR}/src/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/src/tests/PgSAgen-fasta-human.o src/tests/PgSAgen-fasta-human.cpp


${TESTDIR}/src/tests/PgSASearchTest.o: src/tests/PgSASearchTest.cpp 
	${MKDIR} -p ${TESTDIR}/src/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -s -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/src/tests/PgSASearchTest.o src/tests/PgSASearchTest.cpp


${OBJECTDIR}/src/PgSAtest_nomain.o: ${OBJECTDIR}/src/PgSAtest.o src/PgSAtest.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/PgSAtest.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/PgSAtest_nomain.o src/PgSAtest.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/PgSAtest.o ${OBJECTDIR}/src/PgSAtest_nomain.o;\
	fi

${OBJECTDIR}/src/helper_nomain.o: ${OBJECTDIR}/src/helper.o src/helper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/helper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/helper_nomain.o src/helper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/helper.o ${OBJECTDIR}/src/helper_nomain.o;\
	fi

${OBJECTDIR}/src/index/PgSAIndexFactory_nomain.o: ${OBJECTDIR}/src/index/PgSAIndexFactory.o src/index/PgSAIndexFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/index
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/index/PgSAIndexFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/index/PgSAIndexFactory_nomain.o src/index/PgSAIndexFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/index/PgSAIndexFactory.o ${OBJECTDIR}/src/index/PgSAIndexFactory_nomain.o;\
	fi

${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence_nomain.o: ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o src/index/cache/persistence/CountQueriesCachePersistence.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/index/cache/persistence
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence_nomain.o src/index/cache/persistence/CountQueriesCachePersistence.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence.o ${OBJECTDIR}/src/index/cache/persistence/CountQueriesCachePersistence_nomain.o;\
	fi

${OBJECTDIR}/src/pgsaconfig_nomain.o: ${OBJECTDIR}/src/pgsaconfig.o src/pgsaconfig.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pgsaconfig.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pgsaconfig_nomain.o src/pgsaconfig.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pgsaconfig.o ${OBJECTDIR}/src/pgsaconfig_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome_nomain.o: ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o src/pseudogenome/DefaultPseudoGenome.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome_nomain.o src/pseudogenome/DefaultPseudoGenome.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome.o ${OBJECTDIR}/src/pseudogenome/DefaultPseudoGenome_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome_nomain.o: ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o src/pseudogenome/PackedPseudoGenome.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome_nomain.o src/pseudogenome/PackedPseudoGenome.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome.o ${OBJECTDIR}/src/pseudogenome/PackedPseudoGenome_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator_nomain.o: ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator_nomain.o src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.o ${OBJECTDIR}/src/pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator_nomain.o: ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator_nomain.o src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator.o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingDefaultOverlapPseudoGenomeGenerator_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator_nomain.o: ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator_nomain.o src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.o ${OBJECTDIR}/src/pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o: ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility_nomain.o: ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o src/pseudogenome/packing/SymbolsPackingFacility.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/packing
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility_nomain.o src/pseudogenome/packing/SymbolsPackingFacility.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility.o ${OBJECTDIR}/src/pseudogenome/packing/SymbolsPackingFacility_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence_nomain.o: ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o src/pseudogenome/persistence/PseudoGenomePersistence.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/persistence
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence_nomain.o src/pseudogenome/persistence/PseudoGenomePersistence.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence.o ${OBJECTDIR}/src/pseudogenome/persistence/PseudoGenomePersistence_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads_nomain.o: ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o src/pseudogenome/readslist/ListOfConstantLengthReads.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/readslist
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads_nomain.o src/pseudogenome/readslist/ListOfConstantLengthReads.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads.o ${OBJECTDIR}/src/pseudogenome/readslist/ListOfConstantLengthReads_nomain.o;\
	fi

${OBJECTDIR}/src/readsset/DefaultReadsSet_nomain.o: ${OBJECTDIR}/src/readsset/DefaultReadsSet.o src/readsset/DefaultReadsSet.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/readsset
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/readsset/DefaultReadsSet.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/DefaultReadsSet_nomain.o src/readsset/DefaultReadsSet.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/readsset/DefaultReadsSet.o ${OBJECTDIR}/src/readsset/DefaultReadsSet_nomain.o;\
	fi

${OBJECTDIR}/src/readsset/PackedReadsSet_nomain.o: ${OBJECTDIR}/src/readsset/PackedReadsSet.o src/readsset/PackedReadsSet.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/readsset
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/readsset/PackedReadsSet.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/PackedReadsSet_nomain.o src/readsset/PackedReadsSet.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/readsset/PackedReadsSet.o ${OBJECTDIR}/src/readsset/PackedReadsSet_nomain.o;\
	fi

${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator_nomain.o: ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o src/readsset/iterator/ReadsSetIterator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/readsset/iterator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator_nomain.o src/readsset/iterator/ReadsSetIterator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator.o ${OBJECTDIR}/src/readsset/iterator/ReadsSetIterator_nomain.o;\
	fi

${OBJECTDIR}/src/sais/sais_nomain.o: ${OBJECTDIR}/src/sais/sais.o src/sais/sais.c 
	${MKDIR} -p ${OBJECTDIR}/src/sais
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/sais/sais.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -I. -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sais/sais_nomain.o src/sais/sais.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/sais/sais.o ${OBJECTDIR}/src/sais/sais_nomain.o;\
	fi

${OBJECTDIR}/src/suffixarray/DefaultSuffixArray_nomain.o: ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o src/suffixarray/DefaultSuffixArray.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray_nomain.o src/suffixarray/DefaultSuffixArray.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray.o ${OBJECTDIR}/src/suffixarray/DefaultSuffixArray_nomain.o;\
	fi

${OBJECTDIR}/src/suffixarray/SparseSuffixArray_nomain.o: ${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o src/suffixarray/SparseSuffixArray.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/SparseSuffixArray_nomain.o src/suffixarray/SparseSuffixArray.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/suffixarray/SparseSuffixArray.o ${OBJECTDIR}/src/suffixarray/SparseSuffixArray_nomain.o;\
	fi

${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence_nomain.o: ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o src/suffixarray/persistence/SuffixArrayPersistence.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/suffixarray/persistence
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence_nomain.o src/suffixarray/persistence/SuffixArrayPersistence.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence.o ${OBJECTDIR}/src/suffixarray/persistence/SuffixArrayPersistence_nomain.o;\
	fi

${OBJECTDIR}/src/test/PgSAtests_nomain.o: ${OBJECTDIR}/src/test/PgSAtests.o src/test/PgSAtests.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/test
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/test/PgSAtests.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/test/PgSAtests_nomain.o src/test/PgSAtests.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/test/PgSAtests.o ${OBJECTDIR}/src/test/PgSAtests_nomain.o;\
	fi

${OBJECTDIR}/src/test/testdata_nomain.o: ${OBJECTDIR}/src/test/testdata.o src/test/testdata.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/test
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/test/testdata.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -s -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/test/testdata_nomain.o src/test/testdata.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/test/testdata.o ${OBJECTDIR}/src/test/testdata_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f3 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
