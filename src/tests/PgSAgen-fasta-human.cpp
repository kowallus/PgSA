/* 
 * File:   PgSAgen-fasta-human.cpp
 * Author: tomek
 *
 * Created on 2015-02-04, 11:56:48
 */

#include <stdlib.h>
#include <iostream>

#include "../pseudogenome/persistence/PseudoGenomePersistence.h"
#include "../pseudogenome/generator/GreedyVerticalOverlapPseudoGenomeGenerator.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "../readsset/ReadsSetInterface.h"
#include "../readsset/ReadsSetBase.h"
#include "../readsset/DefaultReadsSet.h"
#include "../helper.h"
#include "../suffixarray/DefaultSuffixArrayFactory.h"
#include "../suffixarray/persistence/SuffixArrayPersistence.h"
#include "../suffixarray/generator/SuffixArrayGenerator.h"
#include "../index/PgSAIndexFactory.h"
#include "../index/cache/persistence/CountQueriesCachePersistence.h"

#include "../pseudogenome/generator/PackedPseudoGenomeGenerator.h"
#include "../pseudogenome/packing/SymbolsPackingFacility.h"

//#include "../test/PgSATests.h"

#include <string>

using namespace PgSAIndex;
using namespace PgSAHelpers;
using namespace PgSAReadsSet;

//string readsFile = "2short2"; string idxPrefix = "C://PgSA//2short2";
//string readsFile = "C://PgSA//2sorted_part.txt"; string idxPrefix = "C://PgSA//2sorted_part";
//string readsFile = "C://PgSA//2half.txt"; string idxPrefix = "C://PgSA//2half";
//string readsFile = "//media//pgsadrive//2human"; string idxPrefix = "2human";
string readsFile = "//media//pgsadrive//Human-75-42M.fasta"; string idxPrefix = "3human";

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "PgSAgen-fasta-human test 1" << std::endl;

    PseudoGenomeGeneratorFactory* pggf = new GreedyVerticalOverlapPseudoGenomeGeneratorFactory();

    SuffixArrayBase* sab = SuffixArrayGenerator::generateDefaultSuffixArray(readsFile, pggf);
    
    delete (pggf);
    
     CountQueriesCacheBase* cqcb = sab->getPseudoGenomeBase()->getCountQueriesCacheBase();
    if (cqcb != 0) {
        CountQueriesCachePersistence::writeCache(cqcb, idxPrefix);
        delete (cqcb);   
    }
/**/
    SuffixArrayBase* rsab = SuffixArrayPersistence::readPgSA(idxPrefix + SuffixArrayBase::PGSA_FILE_SUFFIX);
    
    delete (sab);
    delete (rsab);
}

void testSparseGeneration(uchar interval) {
    SuffixArrayBase* rsab = SuffixArrayPersistence::readPgSA(idxPrefix + SuffixArrayBase::PGSA_FILE_SUFFIX);
    SuffixArrayBase* spsab = SuffixArrayGenerator::generateSparseSuffixArray(rsab->getPseudoGenomeBase(), interval);
    SuffixArrayPersistence::writePgSA(spsab, idxPrefix + "_s" + (char) ((int) interval + '0'));
 
    delete(rsab);
    delete(spsab);
}

void testSymbolsPacking() {
    SuffixArrayBase* rsab = SuffixArrayPersistence::readPgSA(idxPrefix + SuffixArrayBase::PGSA_FILE_SUFFIX);
    ReadsSetProperties* prop = rsab->getPseudoGenomeBase()->getReadsSetProperties();

    typedef DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std> DPGType;
    DPGType* dpg = DPGType::castBase(rsab->getPseudoGenomeBase());
    PackedPseudoGenomeGenerator* ppgg = new PackedPseudoGenomeGenerator(dpg, 3);
    typedef PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min> PPGType;
    PPGType* ppg = PPGType::castBase(ppgg->generatePseudoGenomeBase());
    
    cout << dpg->getRead(0) << "\n";
    cout << dpg->getRead(1) << "\n";
    cout << ppg->getRead(0) << "\n";
    cout << ppg->getRead(1) << "\n";
    
    SymbolsPackingFacility<uint_ps_element_min>* spf1 = new SymbolsPackingFacility<uint_ps_element_min>(prop, 2);    
    uchar val = spf1->packSymbols("CT");
    cout << (int) val << ": " << spf1->reverseValue(val,0) << spf1->reverseValue(val,1) << "\n";
    val = spf1->packSymbols("TC");
    cout << (int) val << ": " << spf1->reverseValue(val,0) << spf1->reverseValue(val,1) << "\n";
    val = spf1->clearSuffix(val, 1);
    cout << (int) val << ": " << spf1->reverseValue(val,0) << spf1->reverseValue(val,1) << "\n";

    SymbolsPackingFacility<uint_ps_element_std>* spf2 = new SymbolsPackingFacility<uint_ps_element_std>(prop, 5);    
    ushort sval = spf2->packSymbols("TTCAG");
    cout << (int) sval << ": " << spf2->reverseValue(sval,0) << spf2->reverseValue(sval,1) << spf2->reverseValue(sval,2) << spf2->reverseValue(sval,3) << spf2->reverseValue(sval,4) << "\n";    
    sval = spf2->packSymbols("TTCAT");
    cout << (int) sval << ": " << spf2->reverseValue(sval,0) << spf2->reverseValue(sval,1) << spf2->reverseValue(sval,2) << spf2->reverseValue(sval,3) << spf2->reverseValue(sval,4) << "\n";
    sval = spf2->clearSuffix(sval, 1);
    cout << (int) sval << ": " << spf2->reverseValue(sval,0) << spf2->reverseValue(sval,1) << spf2->reverseValue(sval,2) << spf2->reverseValue(sval,3) << spf2->reverseValue(sval,4) << "\n";
    sval = spf2->packSymbols("GTCAT");
    cout << (int) sval << ": " << spf2->reverseValue(sval,0) << spf2->reverseValue(sval,1) << spf2->reverseValue(sval,2) << spf2->reverseValue(sval,3) << spf2->reverseValue(sval,4) << "\n";
    sval = spf2->clearSuffix(sval, 3);
    cout << (int) sval << ": " << spf2->reverseValue(sval,0) << spf2->reverseValue(sval,1) << spf2->reverseValue(sval,2) << spf2->reverseValue(sval,3) << spf2->reverseValue(sval,4) << "\n";
 
    delete(spf2);
    delete(spf1);
    delete(rsab);
}

void testBoosted() {
    SuffixArrayBase* rsab = SuffixArrayPersistence::readPgSA(idxPrefix + SuffixArrayBase::PGSA_FILE_SUFFIX);

    typedef DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std> DPGType;
    DPGType* dpg = DPGType::castBase(rsab->getPseudoGenomeBase());
    
    cout << dpg->getRead(0) << "\n";
    cout << dpg->getRead(1) << "\n";
    cout << dpg->getRead(2) << "\n";
    cout << dpg->getRead(3) << "\n";
    cout << dpg->getRead(4) << "\n";
    
}
/*
void testIndex() {

    PgSAIndexStandard* pgsa1 = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix);
   
//    vector<uint_reads_cnt_std> r = pgsa1->reportReads("TCGGCAAGGAA");
//    for(int i = 0; i < r.size(); i++) {
//        cout << i << ".\t" << r[i] << "\n";
//    }
    
    string testKmer = "AACGGGACCCT";
    PgSAIndexStandard* pgsa2 = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix + "_s2", idxPrefix); 
    int x = pgsa2->countOccurrences(testKmer);
    cout << x << "\n";
    
    vector<StandardOccurrence> res = pgsa1->reportOccurrences(testKmer);
    for(uint_max i = 0; i < res.size(); i++) {
        cout << i << ".\t" << res[i].first << "\t" << res[i].second << "\n";
    }
    
    
    
    res = pgsa2->reportOccurrences(testKmer);
    for(uint_max i = 0; i < res.size(); i++) {
        cout << i << ".\t" << res[i].first << "\t" << res[i].second << "\n";
    }
    
    
//    StandardOccurrence* a = &res[0];

     cout << pgsa1->reportReads("TCGGCAAGGAA").size() << "\n";
    
    cout << pgsa1->countReads("AATA") << " == " << pgsa1->reportReads("AATA").size() << "\n";
    cout << pgsa1->countReads("AATAA") << " == " << pgsa1->reportReads("AATAA").size() << "\n";
    cout << pgsa1->countReads("AATAATG") << " == " << pgsa1->reportReads("AATAATG").size() << "\n"; 
    cout << pgsa1->countReads("AATAATGCC") << " == " << pgsa1->reportReads("AATAATGCC").size() << "\n"; 
    cout << pgsa1->countReads("AATAATGCCC") << " == " << pgsa1->reportReads("AATAATGCCC").size() << "\n";
    cout << pgsa1->countReads("AATAATGCCCG") << " == " << pgsa1->reportReads("AATAATGCCCG").size() << "\n";

    cout << pgsa1->countOccurrences("AATA") << " == " << pgsa1->reportOccurrences("AATA").size() << "\n";
    cout << pgsa1->countOccurrences("AATAA") << " == " << pgsa1->reportOccurrences("AATAA").size() << "\n";
    cout << pgsa1->countOccurrences("AATAATG") << " == " << pgsa1->reportOccurrences("AATAATG").size() << "\n"; 
    cout << pgsa1->countOccurrences("AATAATGCC") << " == " << pgsa1->reportOccurrences("AATAATGCC").size() << "\n"; 
    cout << pgsa1->countOccurrences("AATAATGCCC") << " == " << pgsa1->reportOccurrences("AATAATGCCC").size() << "\n";
    cout << pgsa1->countOccurrences("AATAATGCCCG") << " == " << pgsa1->reportOccurrences("AATAATGCCCG").size() << "\n";
    
}
*/
void runTestSuite(short interval) {

    string idxFile = idxPrefix;
    if (interval > 1)
        idxFile = idxFile + "_s" + (char) ((int) interval + '0');
    idxFile = idxFile + SuffixArrayBase::PGSA_FILE_SUFFIX;
    
    string cacheFile = idxPrefix + CountQueriesCacheBase::CACHE_FILE_SUFFIX;
    
    PgSAIndexStandard* pgsa1;
    if (interval == 0)
        pgsa1 = PgSAIndexFactory::getPgSAIndexStandard(idxFile, cacheFile, true, true, false);
    else if (interval == -1)
        pgsa1 = PgSAIndexFactory::getPgSAIndexStandard(idxFile, cacheFile, false, false, false);
    else    
        pgsa1 = PgSAIndexFactory::getPgSAIndexStandard(idxFile, cacheFile, true, false, false);

/*    uint_max half = pgsa1->readsCountVirtual() / 2;    
    uint_max i = 0;
    uint_max tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
    i = 1;
    tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
    i = 2;
    tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
    i = 3;
    tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
    i = half;
    tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
    i = pgsa1->readsCountVirtual() - 2;
    tagIdx = (i % 2) * half + i / 2;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";    
    i = pgsa1->readsCountVirtual() - 1;
    tagIdx = (i % 2) * half + i / half;
    cout << i << ". " << pgsa1->getReadVirtual(tagIdx) << "\n";
  */  
//    PgSATest::runTests(pgsa1);
//    PgSATest::runTest(pgsa1, 11, 100000, 22, true);
//    PgSATest::runTest(pgsa1, 11, 100000, 22, false);
    
    delete(pgsa1);

}

/*int main() {

    
    
/*    int bytes = bytesPerSAElement(2400000, 151);
    cout << bytes << "(" << (int) bytesPerValue(2400000) << "," << (int) bytesPerValue(151) << "), " << isSABytesStd(100000000, bytes) << "\n";
*/
/*    
    clock_checkpoint();

    testAutoGeneration();
//    testSparseGeneration(2);    
//    testSparseGeneration(3);
//    testSparseGeneration(4);
//    testSparseGeneration(5);
//    testSparseGeneration(6);
//    testSparseGeneration(7);
//    testSparseGeneration(8);    
   
//    runTestSuite(-1);    
//    runTestSuite(0);
//    runTestSuite(1);
//    runTestSuite(2);      
//    runTestSuite(3);
//    runTestSuite(4);
//    runTestSuite(5);
//    runTestSuite(6);
//    runTestSuite(7);
//    runTestSuite(8);
    
//      testBoosted();
//    testSymbolsPacking();
//    testIndex();
    
    cout << "Time " << clock_millis() << " msec!\n";

}*/

void test2() {
    std::cout << "PgSAgen-fasta-human test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (PgSAgen-fasta-human) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% PgSAgen-fasta-human" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (PgSAgen-fasta-human)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (PgSAgen-fasta-human)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (PgSAgen-fasta-human)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (PgSAgen-fasta-human)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

