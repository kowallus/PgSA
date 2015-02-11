/* 
 * File:   testdata.h
 * Author: Tomek
 *
 * Created on 6 luty 2014, 18:10
 */

#ifndef TESTDATA_H
#define	TESTDATA_H

#include "../pgsaconfig.h"
#include "../helper.h"
#include "../readsset/DefaultReadsSet.h"
#include <random>

using namespace PgSAReadsSet;

namespace PgSATest {

    typedef uint_reads_cnt_std t_reads_c;
    typedef uint_read_len_max t_read;

    std::mt19937 randgenerator;

    std::pair<t_reads_c, t_read>* generateTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount) {

        randgenerator.seed();
        
        std::pair<t_reads_c, t_read>* factors = (std::pair<t_reads_c, t_read>*) malloc(factorsCount * sizeof(std::pair<t_reads_c, t_read>));

        for(int i = 0; i < factorsCount; i++) {
            t_reads_c tagIdx = randgenerator() % readsCount;
            t_read pos = randgenerator() % (readLength - k + 1);
            factors[i].first = tagIdx;
            factors[i].second = pos;
        }

        return factors;
    }

    std::pair<t_reads_c, t_read>* generateSampleTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount) {

        std::pair<t_reads_c, t_read>* factors = (std::pair<t_reads_c, t_read>*) malloc(factorsCount * sizeof(std::pair<t_reads_c, t_read>));

        int pos = 0;
        for(int i = 0; i < factorsCount; i++) {
            factors[i].first = ((pos / (readLength - k + 1)) % readsCount) + 1;
            factors[i].second = pos % (readLength - k + 1);
            pos += (readLength - k + 1) / 2;
        }

        return factors;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    string* generateTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, int factorsCount) {

        string* tFactors = new string[factorsCount];

        std::pair<t_reads_c, t_read>* factors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, factorsCount);

        for(int i = 0; i < factorsCount; i++) {
            t_reads_c tagIdx = factors[i].first;
            t_read pos = factors[i].second;
            tFactors[i] = string(readsSet->getReadVirtual(tagIdx), pos, k);
        }

        free(factors);

        return tFactors;

    }

    template<typename uint_read_len, typename uint_reads_cnt>
    string* generateTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, std::pair<t_reads_c, t_read>* factors, int factorsCount) {

        string* tFactors = new string[factorsCount];

        for(int i = 0; i < factorsCount; i++) {
            t_reads_c tagIdx = factors[i].first;
            t_read pos = factors[i].second;
            tFactors[i] = string(readsSet->getReadVirtual(tagIdx), pos, k);
        }

        return tFactors;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    string* generatePairScrambledTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, int factorsCount) {

        string* tFactors = new string[factorsCount];

        std::pair<t_reads_c, t_read>* factors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, factorsCount);

        t_reads_c half = readsSet->readsCountVirtual() / 2;
        
        for(int i = 0; i < factorsCount; i++) {
            t_reads_c tagIdx = (factors[i].first % 2) * half + factors[i].first / 2;
            t_read pos = factors[i].second;
            tFactors[i] = string(readsSet->getReadVirtual(tagIdx), pos, k);
        }

        free(factors);

        return tFactors;

    }
    
    std::pair<t_reads_c, t_read>* generatePairScrambledTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount) {
        
        std::pair<t_reads_c, t_read>* factors = generateTagFactors(readsCount, readLength, k, factorsCount);
        
        t_reads_c half = readsCount / 2;
        
        for(int i = 0; i < factorsCount; i++) {
            t_reads_c tagIdx = randgenerator() % readsCount;
            t_read pos = randgenerator() % (readLength - k + 1);
            factors[i].first = (tagIdx % 2) * half + tagIdx / 2;
            factors[i].second = pos;
        }

        return factors;
    }
    
}
    
#endif	/* TESTDATA_H */

