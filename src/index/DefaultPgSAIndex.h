#ifndef DEFAULTPGSAINDEX_H_INCLUDED
#define DEFAULTPGSAINDEX_H_INCLUDED

#include "PgSAIndexInterface.h"
#include "../suffixarray/SuffixArrayInterface.h"
#include "cache/DefaultCountQueriesCache.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template <  typename api_uint_reads_cnt, typename api_uint_read_len, // api types
                typename uint_read_len, // types of underlying SA
                typename uint_reads_cnt,
                typename uint_pg_len,
                class SuffixArrayClass>
    class DefaultPgSAIndex: public PgSAIndexInterface< api_uint_reads_cnt, api_uint_read_len >
    {
        private:
            typedef SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass> SuffixArray;
            typedef typename SuffixArray::OccurrencesIterator OccurrencesIterator;
            
            // actually iterator seems to be sufficient
            SuffixArray& sa; 
            
            // cache
            DefaultCountQueriesCacheStandard* countQueriesCache = 0;
            
            string name = "UNKNOWN";
            
        public:

            DefaultPgSAIndex(SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>& sa)
            : sa(sa) {
            };

            DefaultPgSAIndex(SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>& sa,
                    DefaultCountQueriesCacheStandard* countQueriesCache)
            : sa(sa),
              countQueriesCache(countQueriesCache) {
            };
            
            DefaultPgSAIndex(SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>& sa,
                    DefaultCountQueriesCacheStandard* countQueriesCache, const string name)
            : sa(sa),
              countQueriesCache(countQueriesCache),
              name(name) {
            };
            
            virtual ~DefaultPgSAIndex() {
                delete(&sa);
                if (countQueriesCache != 0)
                    delete(countQueriesCache);
            };
            
            // Q1 - In which reads does kmer occur?
            inline vector<api_uint_reads_cnt> reportReadsImpl(const string& kmer) {
                vector<api_uint_reads_cnt> readsIdx;
                reportReadsImpl(kmer, readsIdx);
                return readsIdx;
            };

            inline void reportReadsImpl(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) {
                readsIdxs.clear();
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext()) {
                    if (oit.hasDuplicateFilterFlag()) {
                        if (!oit.hasOccurFlag()){
                            oit.setOccurFlag();
                            readsIdxs.push_back(oit.getReadOriginalIndex());
                        }
                    } else 
                        readsIdxs.push_back(oit.getReadOriginalIndex());
                }
                oit.clearAllOccurFlags();
            };

            inline void reportReadsImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) {
                reportReadsImpl(getKmerByPosition(orgIdx, pos, k), readsIdxs);
            };
            
            // Q2 - In how many reads does kmer occur?
            inline api_uint_reads_cnt countReadsImpl(const string& kmer) {
                if (countQueriesCache != 0) {
                    api_uint_reads_cnt res = countQueriesCache->countReads(kmer);
                    if (res != unknownCacheCount.readsCount)
                        return res;
                }
                
                api_uint_reads_cnt readsCount = 0;
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext()) {
                    if (oit.hasDuplicateFilterFlag()) {
                        if (!oit.hasOccurFlag())
                            oit.setOccurFlag();
                    } else
                        readsCount++;
                }
                return readsCount + oit.clearAllOccurFlags();
            };

            inline api_uint_reads_cnt countReadsImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
                return countReadsImpl(getKmerByPosition(orgIdx, pos, k));
            };

            // Q3 - What are the occurrence positions of kmer?
            inline vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportOccurrencesImpl(const string& kmer) {
                vector<pair<api_uint_reads_cnt, api_uint_read_len>> kmersPos;
                reportOccurrencesImpl(kmer, kmersPos);
                return kmersPos;
            };

            inline void reportOccurrencesImpl(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) {
                kmersPos.clear();
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext())
                    kmersPos.push_back({ oit.getReadOriginalIndex(), oit.getOccurrenceOffset() });
            };

            inline vector<uint_flatten_occurrence_max> reportOccurrencesFlattenImpl(const string& kmer) {
                vector<uint_flatten_occurrence_max> kmersPos;
                reportOccurrencesFlattenImpl(kmer, kmersPos);
                return kmersPos;
            };
            
            inline void reportOccurrencesFlattenImpl(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) {
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext())
                    kmersPos.push_back( oit.getFlattenOccurrence() );
            };
            
            inline void reportOccurrencesImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) {
                reportOccurrencesImpl(getKmerByPosition(orgIdx, pos, k), kmersPos);
            };

            // Q4 - What is the number of occurrences of kmer?
            inline uint_max countOccurrencesImpl(const string& kmer) {
                if (countQueriesCache != 0) {
                    uint_max res = countQueriesCache->countOccurrences(kmer);
                    if (res != unknownCacheCount.occurrencesCount)
                        return res;
                }
                
                uint_max count = 0;
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext())
                    count++;
                
                return count;
            };
            
            inline uint_max countOccurrencesImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
                return countOccurrencesImpl(getKmerByPosition(orgIdx, pos, k));
            };
            
            // Q5 - In which reads does kmer occur only once?
            vector<api_uint_reads_cnt> reportReadsWithSingleOccurrenceImpl(const string& kmer) { 
                vector<api_uint_reads_cnt> readsIdx;
                reportReadsWithSingleOccurrenceImpl(kmer, readsIdx);
                return readsIdx;
            };
            
            void reportReadsWithSingleOccurrenceImpl(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) { 
                readsIdxs.clear();
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext()) {
                    if (oit.hasDuplicateFilterFlag())
                        oit.setOccurOnceFlagAndPushRead();
                    else 
                        readsIdxs.push_back(oit.getReadOriginalIndex());
                }
                
                oit.template clearAllOccurFlagsAndPushReadsWithSingleOccurrence<api_uint_reads_cnt>(readsIdxs);
            };
            
            void reportReadsWithSingleOccurrenceImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) { 
                return reportReadsWithSingleOccurrenceImpl(getKmerByPosition(orgIdx, pos, k), readsIdxs);
            };
            
            // Q6 - In how many reads does kmer occur only once?
            api_uint_reads_cnt countSingleOccurrencesImpl(const string& kmer) { 
                if (countQueriesCache != 0) {
                    api_uint_reads_cnt res = countQueriesCache->countSingleOccurrences(kmer);
                    if (res != unknownCacheCount.singleOccurrencesCount)
                        return res;
                }
                
                api_uint_reads_cnt readsCount = 0;
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext()) {
                    if (oit.hasDuplicateFilterFlag())
                        oit.setOccurOnceFlagAndPushRead();
                    else
                        readsCount++;
                }
                return readsCount + oit.clearAllOccurFlagsAndCountSingleOccurrences();
            };
            
            api_uint_reads_cnt countSingleOccurrencesImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) { 
                return countSingleOccurrencesImpl(getKmerByPosition(orgIdx, pos, k));
            };
            
            // Q7 - What are the occurrence positions of kmer in the reads where it occurs only once?
            vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportSingleOccurrencesImpl(const string& kmer) { 
                vector<pair<api_uint_reads_cnt, api_uint_read_len>> kmersPos;
                reportSingleOccurrencesImpl(kmer, kmersPos);
                return kmersPos;
            };
            
            void reportSingleOccurrencesImpl(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { 
                kmersPos.clear();
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext())
                    if (oit.hasDuplicateFilterFlag())
                        oit.setOccurOnceFlagAndPushOccurrence();
                    else 
                        kmersPos.push_back({ oit.getReadOriginalIndex(), oit.getOccurrenceOffset() });
                
                oit.template clearAllOccurFlagsAndPushSingleOccurrences<api_uint_reads_cnt, api_uint_read_len>(kmersPos);
            };
            
            vector<uint_flatten_occurrence_max> reportSingleOccurrencesFlattenImpl(const string& kmer) { 
                vector<uint_flatten_occurrence_max> kmersPos;
                reportSingleOccurrencesFlattenImpl(kmer, kmersPos);
                return kmersPos;
            };
            
            void reportSingleOccurrencesFlattenImpl(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) { 
                kmersPos.clear();
                OccurrencesIterator& oit = sa.getKmerOccurrencesIterator(kmer);
                while(oit.moveNext())
                    if (oit.hasDuplicateFilterFlag())
                        oit.setOccurOnceFlagAndPushOccurrence();
                    else 
                        kmersPos.push_back( oit.getFlattenOccurrence() );
                
                oit.clearAllOccurFlagsAndPushSingleOccurrencesFlatten(kmersPos);
            };
            
            void reportSingleOccurrencesImpl(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { 
                return reportSingleOccurrencesImpl(getKmerByPosition(orgIdx, pos, k), kmersPos);
            };
           
            const string getDescription() {
                string desc = name + "\n" + sa.getDescription();
                if (countQueriesCache != 0)
                    desc = desc + countQueriesCache->getDescription();
                
                return desc;
            }
            
            inline string getKmerByPosition(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
                return sa.getKmer(orgIdx, pos, k);
            }
            
            bool isReadLengthConstantVirtual() {
                return sa.getReadsSet()->isReadLengthConstantVirtual();
            };
            
            api_uint_read_len maxReadLengthVirtual() {
                return sa.getReadsSet()->maxReadLengthVirtual();
            };

            api_uint_reads_cnt readsCountVirtual() {
                return sa.getReadsSet()->readsCountVirtual();
            };

            const string getReadVirtual(api_uint_reads_cnt originalIdx) {
                return sa.getReadsSet()->getReadVirtual(originalIdx);
            };
            
            api_uint_read_len readLengthVirtual(api_uint_reads_cnt originalIdx) {
                return sa.getReadsSet()->readLengthVirtual(originalIdx);
            };
            
            // Q1 - In which reads does kmer occur?
            vector<api_uint_reads_cnt> reportReads(const string& kmer) { return reportReadsImpl(kmer); };
            void reportReads(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) { reportReadsImpl(kmer, readsIdxs); };
            void reportReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) { reportReadsImpl(orgIdx, pos, k, readsIdxs); };
            // Q2 - In how many reads does kmer occur?
            api_uint_reads_cnt countReads(const string& kmer) { return countReadsImpl(kmer); };
            api_uint_reads_cnt countReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) { return countReadsImpl(orgIdx, pos, k); };
            // Q3 - What are the occurrence positions of kmer?
            vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportOccurrences(const string& kmer) { return reportOccurrencesImpl(kmer); };
            void reportOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& readsIdxs) { reportOccurrencesImpl(kmer, readsIdxs); };
            vector<uint_flatten_occurrence_max> reportOccurrencesFlatten(const string& kmer) { return reportOccurrencesFlattenImpl(kmer); };
            void reportOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& readsIdxs) { reportOccurrencesFlattenImpl(kmer, readsIdxs); }
            void reportOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& readsIdxs) { reportOccurrencesImpl(orgIdx, pos, k, readsIdxs); };
            // Q4 - What is the number of occurrences of kmer?
            uint_max countOccurrences(const string& kmer) { return countOccurrencesImpl(kmer); };
            uint_max countOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) { return countOccurrencesImpl(orgIdx, pos, k); };          
            // Q5 - In which reads does kmer occur only once?
            vector<api_uint_reads_cnt> reportReadsWithSingleOccurrence(const string& kmer) { return reportReadsWithSingleOccurrenceImpl(kmer); };
            void reportReadsWithSingleOccurrence(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) { reportReadsWithSingleOccurrenceImpl(kmer, readsIdxs); };
            void reportReadsWithSingleOccurrence(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) { reportReadsWithSingleOccurrenceImpl(orgIdx, pos, k, readsIdxs); };           
            // Q6 - In how many reads does kmer occur only once?
            api_uint_reads_cnt countSingleOccurrences(const string& kmer) { return countSingleOccurrencesImpl(kmer); };
            api_uint_reads_cnt countSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) { return countSingleOccurrencesImpl(orgIdx, pos, k); };            
            // Q7 - What are the occurrence positions of kmer in the reads where it occurs only once?
            vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportSingleOccurrences(const string& kmer) { return reportSingleOccurrencesImpl(kmer); };
            void reportSingleOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { reportSingleOccurrencesImpl(kmer, kmersPos); };
            vector<uint_flatten_occurrence_max> reportSingleOccurrencesFlatten(const string& kmer) { return reportSingleOccurrencesFlattenImpl(kmer); };
            void reportSingleOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) { reportSingleOccurrencesFlattenImpl(kmer, kmersPos); };
            void reportSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { reportSingleOccurrencesImpl(orgIdx, pos, k, kmersPos); };            
            
    };

}

#endif // DEFAULTPGSAINDEX_H_INCLUDED
