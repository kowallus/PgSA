#ifndef PGSAINDEXINTERFACE_H_INCLUDED
#define PGSAINDEXINTERFACE_H_INCLUDED

#include "../pgsaconfig.h"
#include "../readsset/ReadsSetInterface.h"

#include <vector>

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename api_uint_reads_cnt, typename api_uint_read_len >
    class PgSAIndexInterface: public ReadsSetInterface<api_uint_read_len, api_uint_reads_cnt>
    {
        public:

            virtual ~PgSAIndexInterface() {};
            
            // Q1 - In which reads does kmer occur?
            virtual vector<api_uint_reads_cnt> reportReads(const string& kmer) = 0;
            virtual void reportReads(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) = 0;

            virtual void reportReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) = 0;
            
            // Q2 - In how many reads does kmer occur?
            virtual api_uint_reads_cnt countReads(const string& kmer) = 0;

            virtual api_uint_reads_cnt countReads(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) = 0;
            
            // Q3 - What are the occurrence positions of kmer?
            virtual vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportOccurrences(const string& kmer) = 0;
            virtual void reportOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) = 0;
            virtual vector<uint_flatten_occurrence_max> reportOccurrencesFlatten(const string& kmer) = 0;
            virtual void reportOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) = 0;

            virtual void reportOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) = 0;
            
            // Q4 - What is the number of occurrences of kmer?
            virtual uint_max countOccurrences(const string& kmer) = 0;
            
            virtual uint_max countOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) = 0;
            
            // Q5 - In which reads does kmer occur only once?
            virtual vector<api_uint_reads_cnt> reportReadsWithSingleOccurrence(const string& kmer) = 0;
            virtual void reportReadsWithSingleOccurrence(const string& kmer, vector<api_uint_reads_cnt>& readsIdxs) = 0;
            
            virtual void reportReadsWithSingleOccurrence(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<api_uint_reads_cnt>& readsIdxs) = 0;
            
            // Q6 - In how many reads does kmer occur only once?
            virtual api_uint_reads_cnt countSingleOccurrences(const string& kmer) = 0;
            
            virtual api_uint_reads_cnt countSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k) = 0;
            
            // Q7 - What are the occurrence positions of kmer in the reads where it occurs only once?
            virtual vector<pair<api_uint_reads_cnt, api_uint_read_len>> reportSingleOccurrences(const string& kmer) = 0;
            virtual void reportSingleOccurrences(const string& kmer, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) = 0;
            virtual vector<uint_flatten_occurrence_max> reportSingleOccurrencesFlatten(const string& kmer) = 0;
            virtual void reportSingleOccurrencesFlatten(const string& kmer, vector<uint_flatten_occurrence_max>& kmersPos) = 0;

            virtual void reportSingleOccurrences(api_uint_reads_cnt orgIdx, api_uint_read_len pos, api_uint_read_len k, vector<pair<api_uint_reads_cnt, api_uint_read_len>>& kmersPos) = 0;
            
            virtual const string getDescription() = 0;
            
    };

    typedef PgSAIndexInterface<uint_reads_cnt_std, unsigned int> PgSAIndexStandard;
    typedef ReadsSetInterface<unsigned int, uint_reads_cnt_std> ReadsSetStandard;
    typedef pair< uint_reads_cnt_std, unsigned int > StandardOccurrence;

}

#endif // PGSAINDEXINTERFACE_H_INCLUDED
