/* 
 * File:   OccurrencesIteratorInterface.h
 * Author: Tomek
 *
 * Created on 3 luty 2014, 20:32
 */

#ifndef OCCURRENCESITERATORINTERFACE_H
#define	OCCURRENCESITERATORINTERFACE_H

#include "../../pseudogenome/readslist/iterator/ReadsListIteratorInterface.h"

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt, class ReadsListIteratorClass, class OccurrencesIteratorClass >
    class OccurrencesIteratorInterface: public ReadsIteratorInterface<uint_read_len, uint_reads_cnt, OccurrencesIteratorClass>
    {
        protected:
            typedef ReadsListIteratorInterface<uint_read_len, uint_reads_cnt, ReadsListIteratorClass> ReadsListIterator;
            
            ReadsListIterator& readsIterator;
        
        public:
            
            OccurrencesIteratorInterface(ReadsListIteratorClass& readsIterator)
            : readsIterator(readsIterator)
            {
            };
            
            virtual ~OccurrencesIteratorInterface() {};

            inline void findOccurrencesOf(const char* kmer, const uint_read_len kmerLength) { static_cast<OccurrencesIteratorClass*>(this)->findOccurrencesOfImpl(kmer, kmerLength); };
            
            // delegates
            inline uint_reads_cnt getReadOriginalIndexImpl() { return readsIterator.getReadOriginalIndex(); };
            inline uint_read_len getOccurrenceOffsetImpl() { return readsIterator.getOccurrenceOffset(); };
            inline uint_flatten_occurrence_max getFlattenOccurrenceImpl() { return readsIterator.getFlattenOccurrence(); };
            
            inline bool hasDuplicateFilterFlagImpl() { return readsIterator.hasDuplicateFilterFlag(); };
            inline bool hasOccurFlagImpl() { return readsIterator.hasOccurFlag(); };
            inline bool hasOccurOnceFlagImpl() { return readsIterator.hasOccurOnceFlag(); };
            inline void setOccurFlagImpl() { readsIterator.setOccurFlag(); };            
            
            inline uint_reads_cnt clearAllOccurFlagsImpl() { return readsIterator.clearAllOccurFlags(); };
            
            inline void setOccurOnceFlagAndPushReadImpl() { readsIterator.setOccurOnceFlagAndPushRead(); };
            inline void setOccurOnceFlagAndPushOccurrenceImpl() { readsIterator.setOccurOnceFlagAndPushOccurrence(); };
            
            template<typename api_uint_reads_cnt>
            inline void clearAllOccurFlagsAndPushReadsWithSingleOccurrenceImpl(vector<api_uint_reads_cnt>& reads) { readsIterator.template clearAllOccurFlagsAndPushReadsWithSingleOccurrence<api_uint_reads_cnt>(reads);  };
            
            inline uint_reads_cnt clearAllOccurFlagsAndCountSingleOccurrencesImpl() { return readsIterator.clearAllOccurFlagsAndCountSingleOccurrences(); };
            
            template<typename api_uint_reads_cnt, typename api_uint_read_len>
            inline void clearAllOccurFlagsAndPushSingleOccurrencesImpl(vector<pair<api_uint_reads_cnt, api_uint_read_len>>& occurrences) { readsIterator.template clearAllOccurFlagsAndPushSingleOccurrences<api_uint_reads_cnt, api_uint_read_len>(occurrences);  };
                        
            inline void clearAllOccurFlagsAndPushSingleOccurrencesFlattenImpl(vector<uint_flatten_occurrence_max>& flattenOccurrences) { readsIterator.clearAllOccurFlagsAndPushSingleOccurrencesFlatten(flattenOccurrences);  };            
    };

}


#endif	/* OCCURRENCESITERATORINTERFACE_H */

