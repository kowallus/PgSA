/* 
 * File:   gkbugdump.cpp
 * Author: tomek
 *
 * Created on 2015-02-11, 09:10:19
 */

#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "../index/PgSAIndexFactory.h"
#include "../helper.h"

#include <libGkArrays/gkArrays.h>

/*
 * Simple C++ Test Suite
 */

string idxPrefix = "//media//pgsadrive//PgSA-test//2pseudo";
char* readsFile = "/media/pgsadrive/Ecoli_R.fastq";
int readLength = 151;


void testGkBug() {
    std::cout << "gkbugdump GkBug test 1" << std::endl;
    PgSAIndexStandard* pgsaIndex = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix + ".pgsa", idxPrefix + ".pgc", false);

    vector<StandardOccurrence> q3res_pgsa;
    pgsaIndex->reportOccurrences("AGTTGAACTGC", q3res_pgsa);
    std::cout << "|Q3(AGTTGAACTGC)| = " << q3res_pgsa.size() << std::endl;
    
    delete pgsaIndex;

    int k = 11;
    
    gkarrays::gkArrays *reads;

    std::pair <uint, uint> * q3res_gka;
    uint nb_fact;
    
    reads = new gkarrays::gkArrays(readsFile, k, false, readLength, true);

    q3res_gka = reads->getTagsWithFactor("AGTTGAACTGC", k, nb_fact);
    
    std::cout << "|Q3(AGTTGAACTGC)| = " << nb_fact << std::endl;
    
    delete reads;

    std::vector<uint> pos_pgsa;
    std::vector<uint> pos_gk;
    
    for (std::vector<StandardOccurrence>::iterator it = q3res_pgsa.begin() ; it != q3res_pgsa.end(); ++it)
        pos_pgsa.push_back(it->first * readLength + it->second); 
    std::sort(pos_pgsa.begin(), pos_pgsa.end());
   
    for (int i = 0; i < nb_fact; ++i)
        pos_gk.push_back(q3res_gka[i].first * readLength + q3res_gka[i].second);
    std::sort(pos_gk.begin(), pos_gk.end());
        
    std::vector<uint>::iterator it_pgsa;
    std::vector<uint>::iterator it_gk;
    for (it_pgsa = pos_pgsa.begin(), it_gk = pos_gk.begin(); 
            it_pgsa != pos_pgsa.end() && it_gk != pos_gk.end(); ++it_pgsa, ++it_gk)
        std::cout << "PgSA/GkA:\t" << *it_pgsa << "\t/\t" << *it_gk << std::endl;
    
    for (;it_pgsa != pos_pgsa.end(); ++it_pgsa)
        std::cout << "PgSA:\t" << *it_pgsa << std::endl;
    
    
    for (;it_gk != pos_gk.end(); ++it_gk)
        std::cout << "GkA:\t" << *it_gk << std::endl;
}

int k = 11;

bool scan3797557(gkarrays::gkArrays *reads, char* factor) {
    
    std::pair <uint, uint> * q3res_gka; 
    uint nb_fact;    

    q3res_gka = reads->getTagsWithFactor(factor, k, nb_fact);
    
    std::cout << "|Q3(" << string(factor) << ")| = " << nb_fact << std::endl;
    
    for (int i = 0; i < nb_fact; ++i)
        if (37977557 == q3res_gka[i].first * readLength + q3res_gka[i].second)
            cout << "Q3(" << string(factor) << ") contains 37977557" << std::endl;
    
}

void testGkNk() {
    std::cout << "gkbugdump GkNk test 1" << std::endl;
    PgSAIndexStandard* pgsaIndex = PgSAIndexFactory::getPgSAIndexStandard(idxPrefix + ".pgsa", idxPrefix + ".pgc", false);

    std::cout << "|Q3(AGTTGAACTGC)| = " << pgsaIndex->countOccurrences("AGTTGAACTGC") << std::endl;
    std::cout << "|Q3(CGTTGAACTGC)| = " << pgsaIndex->countOccurrences("CGTTGAACTGC") << std::endl;
    std::cout << "|Q3(GGTTGAACTGC)| = " << pgsaIndex->countOccurrences("GGTTGAACTGC") << std::endl;
    std::cout << "|Q3(TGTTGAACTGC)| = " << pgsaIndex->countOccurrences("TGTTGAACTGC") << std::endl;
    std::cout << "|Q3(NGTTGAACTGC)| = " << pgsaIndex->countOccurrences("NGTTGAACTGC") << std::endl;
    
    delete pgsaIndex;

    gkarrays::gkArrays *reads;
    
    reads = new gkarrays::gkArrays(readsFile, k, false, readLength, true);

    scan3797557(reads, "AGTTGAACTGC");
    scan3797557(reads, "CGTTGAACTGC");
    scan3797557(reads, "GGTTGAACTGC");
    scan3797557(reads, "TGTTGAACTGC");
    scan3797557(reads, "NGTTGAACTGC");
    
    delete reads;
    
  
}

int main(int argc, char** argv) {
    
//    testGkBug();
    
    testGkNk();
    
    return (EXIT_SUCCESS);
}

