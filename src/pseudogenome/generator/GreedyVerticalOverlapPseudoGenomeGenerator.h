#ifndef GREEDYVERTICALOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define GREEDYVERTICALOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "../../index/cache/DefaultCountQueriesCache.h"
#include "AbstractOverlapPseudoGenomeGenerator.h"
#include <algorithm>
#include <set>

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt >
    class GreedyVerticalOverlapGeneratorTemplate: public AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>
    {
        private:

            // auxiliary structures
            uint_reads_cnt dupsTotal;

            DefaultReadsSet* orgReadsSet = 0;
            
            virtual uint_read_len readLength(uint_reads_cnt incIdx) override;
            virtual string getReadUpToOverlap(uint_reads_cnt incIdx) override;
            virtual uint_reads_cnt readsTotal() override;
            
            virtual ReadsSetProperties* getReadsSetProperties() override;
            
            struct ReadWithIdx{
                uint_reads_cnt incIdx;
                const char* read;
            };

            struct ReadsSufPreComparator {
                bool operator() (ReadWithIdx suffix, ReadWithIdx prefix) const {
                    return readsSufPreCmp(suffix.read, prefix.read) < 0;
                }
            };

            std::multiset<ReadWithIdx, ReadsSufPreComparator> readsSet;

            void populateReadsSet();
            uint_reads_cnt matchInReadsSet(uint_reads_cnt suffixIdx, uint_read_len overlap);
            void overlapMultisetVertical(uint_read_len, uint_read_len);

            virtual void findOverlappingReads() override;
            
        public:

            GreedyVerticalOverlapGeneratorTemplate(DefaultReadsSet* readsSet);
            virtual ~GreedyVerticalOverlapGeneratorTemplate();

            bool isPseudoGenomeLengthStandardVirtual();
            bool isPseudoGenomeLengthMaximalVirtual();

    };

    class GreedyVerticalOverlapPseudoGenomeGeneratorFactory: public PseudoGenomeGeneratorFactory
    {
        private:

            template<typename uint_read_len, typename uint_reads_cnt>
            PseudoGenomeGeneratorBase* getGeneratorFullTemplate(DefaultReadsSet* readsSet);

            template<typename uint_read_len>
            PseudoGenomeGeneratorBase* getGeneratorPartialTemplate(DefaultReadsSet* readsSet);

        public:

            GreedyVerticalOverlapPseudoGenomeGeneratorFactory() {};

            PseudoGenomeGeneratorBase* getGenerator(string readsFile, string pairFile);

    };

}


#endif // GREEDYVERTICALOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
