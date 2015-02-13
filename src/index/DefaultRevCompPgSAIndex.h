/* 
 * File:   DefaultRevCompPgSAIndex.h
 * Author: Tomek
 *
 * Created on 11 marzec 2014, 16:37
 */

#ifndef DEFAULTREVCOMPPGSAINDEX_H
#define	DEFAULTREVCOMPPGSAINDEX_H

#include "DefaultPgSAIndex.h"
#include <exception>

using namespace PgSAHelpers;

namespace PgSAIndex {

    template <  typename api_uint_reads_cnt, typename api_uint_read_len, // api types
                    typename uint_read_len, // types of underlying SA
                    typename uint_reads_cnt,
                    typename uint_pg_len,
                    class SuffixArrayClass>
    class DefaultRevCompPgSAIndex: public PgSAIndexInterface< api_uint_reads_cnt, api_uint_read_len > {
        private:

            typedef DefaultPgSAIndex<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass> DefaultPgSAImpl;


            DefaultPgSAImpl* impl;

        public:

            DefaultRevCompPgSAIndex(DefaultPgSAImpl* impl)
            : impl(impl) {};


            ~DefaultRevCompPgSAIndex() {
                delete(impl);
            };

            // Q1 - In which reads does kmer occur?
            vector<api_uint_reads_cnt> reportReads(const string& kmer) { 
                vector<api_uint_reads_cnt> readsIdx;
                reportReads(kmer, readsIdx);
                return readsIdx;
            };

            void reportReads(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) { 
                impl->reportReadsImpl(kmer, readsIdxs);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportReadsImpl(revComp, readsIdxs);
            };
            
            void reportReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) {
//                throw std::runtime_error("not implemented");
            }

            // Q2 - In how many reads does kmer occur?
            api_uint_reads_cnt countReads(const string& kmer) { 
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    return impl->countReadsImpl(kmer) + impl->countReadsImpl(reverseComplement(kmer));
                else 
                    return impl->countReadsImpl(kmer);
            };

            api_uint_reads_cnt countReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
//                throw std::runtime_error("not implemented");
            }
            
            // Q3 - What are the occurrence positions of kmer?
            vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportOccurrences(const string& kmer) { 
                vector<pair<api_uint_reads_cnt, api_uint_read_len>> kmersPos;
                reportOccurrences(kmer, kmersPos);
                return kmersPos;
            };

            void reportOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { 
                impl->reportOccurrencesImpl(kmer, kmersPos);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportOccurrencesImpl(revComp, kmersPos);
            };

            vector<uint_flatten_occurrence_max> reportOccurrencesFlatten(const string& kmer) { 
                vector<uint_flatten_occurrence_max> kmersPos;
                reportOccurrencesFlatten(kmer, kmersPos);
                return kmersPos;
            };


            void reportOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) { 
                impl->reportOccurrencesFlattenImpl(kmer, kmersPos);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportOccurrencesFlattenImpl(revComp, kmersPos);
            };

            void reportOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) {
//                throw std::runtime_error("not implemented");
            }
            
            // Q4 - What is the number of occurrences of kmer?
            uint_max countOccurrences(const string& kmer) { 
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    return impl->countOccurrences(kmer) + impl->countOccurrences(reverseComplement(kmer));
                else
                    return impl->countOccurrences(kmer);
            };

            uint_max countOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
//                throw std::runtime_error("not implemented");
            }
            
            // Q5 - In which reads does kmer occur only once?
            vector<api_uint_reads_cnt> reportReadsWithSingleOccurrence(const string& kmer) { 
                vector<api_uint_reads_cnt> readsIdx;
                reportReadsWithSingleOccurrence(kmer, readsIdx);
                return readsIdx;
            };

            void reportReadsWithSingleOccurrence(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) { 
                impl->reportReadsWithSingleOccurrenceImpl(kmer, readsIdxs);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportReadsWithSingleOccurrenceImpl(revComp, readsIdxs);
            };

            void reportReadsWithSingleOccurrence(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) {
//                throw std::runtime_error("not implemented");
            }
            
            // Q6 - In how many reads does kmer occur only once?
            api_uint_reads_cnt countSingleOccurrences(const string& kmer) { 
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    return impl->countSingleOccurrences(kmer) + impl->countSingleOccurrences(reverseComplement(kmer));
                else
                    return impl->countSingleOccurrences(kmer);
            };
            
            api_uint_reads_cnt countSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) {
//                throw std::runtime_error("not implemented");
            }

            // Q7 - What are the occurrence positions of kmer in the reads where it occurs only once?
            vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportSingleOccurrences(const string& kmer) { 
                vector<pair<api_uint_reads_cnt, api_uint_read_len>> kmersPos;
                reportSingleOccurrences(kmer, kmersPos);
                return kmersPos;
            };

            void reportSingleOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) { 
                impl->reportSingleOccurrencesImpl(kmer, kmersPos);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportSingleOccurrencesImpl(revComp, kmersPos);
            };

            vector<uint_flatten_occurrence_max> reportSingleOccurrencesFlatten(const string& kmer) { 
                vector<uint_flatten_occurrence_max> kmersPos;
                reportSingleOccurrencesFlatten(kmer, kmersPos);
                return kmersPos;        
            };

            void reportSingleOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) { 
                impl->reportSingleOccurrencesFlattenImpl(kmer, kmersPos);
                string revComp = reverseComplement(kmer);
                if (revComp != kmer)
                    impl->reportSingleOccurrencesFlattenImpl(revComp, kmersPos);
            };

            void reportSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) {
//                throw std::runtime_error("not implemented");
            }
            
            const string getDescription() { return impl->getDescription(); };

            bool isReadLengthConstantVirtual() { return impl->isReadLengthConstantVirtual(); };

            api_uint_read_len maxReadLengthVirtual() { return impl->maxReadLengthVirtual(); };

            api_uint_reads_cnt readsCountVirtual() { return impl->readsCountVirtual(); };

            const string getReadVirtual(api_uint_reads_cnt originalIdx) { return impl->getReadVirtual(originalIdx); };

            api_uint_read_len readLengthVirtual(api_uint_reads_cnt originalIdx) { return impl->readLengthVirtual(originalIdx); };

    };
}
    
#endif	/* DEFAULTREVCOMPPGSAINDEX_H */

