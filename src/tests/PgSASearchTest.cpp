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
string fixKIdxPrefix = "d:/Glony/pgsatest/2pseudoidx";
//string fixKIdxPrefix = "d:/Glony/pgsatest/2pseudoidx_k22";
int readLength = 151;


void testSearch() {
    std::cout << "PgSASearchTest test 1" << std::endl;
    PgSAIndexStandard* pgsaIndex = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix + ".pgsa", "", false);
    PgSAIndexStandard* pgsaIndexFixK = PgSAIndexFactory::getPgSAIndexStandard(fixKIdxPrefix + ".pgsa", "", false);
    
    //string kmer = "TGCACCTGCGGTGGCCGATAAAGCCGAC";
    string kmer = "TTTTTTTTTTTTTTTTTTTTTT";
    //string kmer = "TGCACCTGCGGTGGCCGATAAA";
    //  string kmer = "ATGGCGTTAATCAGCGAAGTTT";
    
    vector<StandardOccurrence> q3res_pgsa;
    vector<StandardOccurrence> q3res_pgsa_fixk;

    pgsaIndex->reportOccurrences(kmer, q3res_pgsa);
    std::cout << "|Q3(" << kmer << ")| = " << q3res_pgsa.size() << std::endl;
    pgsaIndexFixK->reportOccurrences(kmer, q3res_pgsa_fixk);    
    std::cout << "|Q3(" << kmer << ")| = " << q3res_pgsa_fixk.size() << std::endl;
    
    int limit = q3res_pgsa.size();
    if (limit > q3res_pgsa_fixk.size())
        limit = q3res_pgsa_fixk.size();
    int j = 0;
    
    for(int i = 0; i < limit; i++, j++) {
        StandardOccurrence occ1 = q3res_pgsa[i];
        StandardOccurrence occ2 = q3res_pgsa_fixk[j];
//        if(occ1.first != occ2.first || occ1.second != occ2.second) {
            cout << i << " - " << j << ": " << occ1.first << " - " << occ2.first << " and " << occ1.second << " - " << occ2.second << std::endl;
//            j--;
//        }
    }
    StandardOccurrence occLast = q3res_pgsa_fixk[1655];
    cout << occLast.first << " " << occLast.second << std::endl;
    //TGCACCTGCGGTGGCCGATAAAGCCGAC 

    
    delete pgsaIndex;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% PgSASearchTest" << std::endl;
    testSearch();

    return (EXIT_SUCCESS);
}

