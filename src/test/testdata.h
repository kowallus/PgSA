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

using namespace PgSAReadsSet;

namespace PgSATest {

    typedef uint_reads_cnt_std t_reads_c;
    typedef uint_read_len_max t_read;

    
    std::pair<t_reads_c, t_read>* generateTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount);

    std::pair<t_reads_c, t_read>* generateSampleTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount);

    template<typename uint_read_len, typename uint_reads_cnt>
    string* generateTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, int factorsCount);
    template<typename uint_read_len, typename uint_reads_cnt>
    string* generateTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, std::pair<t_reads_c, t_read>* factors, int factorsCount);
    
    template<typename uint_read_len, typename uint_reads_cnt>
    string* generatePairScrambledTextFactors(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, t_read k, int factorsCount);
    
    std::pair<t_reads_c, t_read>* generatePairScrambledTagFactors(t_reads_c readsCount, t_read readLength, t_read k, int factorsCount);     
            
}
    
#endif	/* TESTDATA_H */

