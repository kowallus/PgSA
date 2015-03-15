/* 
 * File:   PgSASearchTest.cpp
 * Author: coach
 *
 * Created on 2015-03-11, 13:51:16
 */

#include <stdlib.h>
#include <iostream>

#include "../index/PgSAIndexFactory.h"
#include "../helper.h"

/*
 * Simple C++ Test Suite
 */

string idxPrefix = "d:/kov/noe2-testy2015/PgSA-s/2pseudo";
int readLength = 151;


void testSearch() {
    std::cout << "PgSASearchTest test 1" << std::endl;
    PgSAIndexStandard* pgsaIndex = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix + ".pgsa", idxPrefix + ".pgc", false);

    vector<StandardOccurrence> q3res_pgsa;
    pgsaIndex->reportOccurrences("TGCACCTGCGGTGGCCGATAAAGCCGAC", q3res_pgsa);
    std::cout << "|Q3(TGCACCTGCGGTGGCCGATAAAGCCGAC)| = " << q3res_pgsa.size() << std::endl;
    
    //TGCACCTGCGGTGGCCGATAAAGCCGAC 

    
    delete pgsaIndex;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% PgSASearchTest" << std::endl;
    testSearch();

    return (EXIT_SUCCESS);
}

