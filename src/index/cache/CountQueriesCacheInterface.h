/* 
 * File:   CountQueriesCacheInterface.h
 * Author: Tomek
 *
 * Created on 11 luty 2014, 14:53
 */

#include "CountQueriesCacheBase.h"

#ifndef COUNTQUERIESCACHEINTERFACE_H
#define	COUNTQUERIESCACHEINTERFACE_H

namespace PgSAIndex {

    template<typename api_uint_reads_cnt, class CountQueriesCacheClass>
    class CountQueriesCacheInterface: public CountQueriesCacheBase
    {
        public:

            // Q2 - In how many reads does kmer occur?
            inline api_uint_reads_cnt countReads(const string& kmer) { return static_cast<CountQueriesCacheClass*>(this)->countReadsImpl(kmer); };

            // Q4 - What is the number of occurrences of kmer?
            inline uint_max countOccurrences(const string& kmer) { return static_cast<CountQueriesCacheClass*>(this)->countOccurrencesImpl(kmer); };

            // Q6 - In how many reads does kmer occur only once?
            inline api_uint_reads_cnt countSingleOccurrences(const string& kmer) { return static_cast<CountQueriesCacheClass*>(this)->countSingleOccurrencesImpl(kmer); };

            inline const string getDescription() { return static_cast<CountQueriesCacheClass*>(this)->getDescriptionImpl(); };
            
            static CountQueriesCacheInterface<api_uint_reads_cnt, CountQueriesCacheClass>* castBase(CountQueriesCacheBase* base) {
                // TODO: validate
                return static_cast<CountQueriesCacheInterface<api_uint_reads_cnt, CountQueriesCacheClass>*>(base);
            }
    };

}

#endif	/* COUNTQUERIESCACHEINTERFACE_H */

