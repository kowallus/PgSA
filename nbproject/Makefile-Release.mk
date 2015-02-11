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
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include MakefileNetBeans

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/helper.o \
	${OBJECTDIR}/src/pgsaconfig.o \
	${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f2

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pipe
CXXFLAGS=-pipe

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libPgSA.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libPgSA.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libPgSA.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/helper.o: src/helper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/helper.o src/helper.cpp

${OBJECTDIR}/src/pgsaconfig.o: src/pgsaconfig.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pgsaconfig.o src/pgsaconfig.cpp

${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o: src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/src/tests/PgSAgen-fasta-human.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/f2: ${TESTDIR}/src/tests/gkbugdump.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} 


${TESTDIR}/src/tests/PgSAgen-fasta-human.o: src/tests/PgSAgen-fasta-human.cpp 
	${MKDIR} -p ${TESTDIR}/src/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/src/tests/PgSAgen-fasta-human.o src/tests/PgSAgen-fasta-human.cpp


${TESTDIR}/src/tests/gkbugdump.o: src/tests/gkbugdump.cpp 
	${MKDIR} -p ${TESTDIR}/src/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Wall -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/src/tests/gkbugdump.o src/tests/gkbugdump.cpp


${OBJECTDIR}/src/helper_nomain.o: ${OBJECTDIR}/src/helper.o src/helper.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/helper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/helper_nomain.o src/helper.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/helper.o ${OBJECTDIR}/src/helper_nomain.o;\
	fi

${OBJECTDIR}/src/pgsaconfig_nomain.o: ${OBJECTDIR}/src/pgsaconfig.o src/pgsaconfig.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pgsaconfig.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pgsaconfig_nomain.o src/pgsaconfig.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pgsaconfig.o ${OBJECTDIR}/src/pgsaconfig_nomain.o;\
	fi

${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o: ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/pseudogenome/generator
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O3 -Wall -std=c++11 -fPIC  -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.o ${OBJECTDIR}/src/pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libPgSA.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
