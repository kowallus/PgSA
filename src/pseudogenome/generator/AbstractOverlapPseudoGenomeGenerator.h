#ifndef ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "../../index/cache/DefaultCountQueriesCache.h"
#include <algorithm>
#include <set>

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt >
    class AbstractOverlapPseudoGenomeGeneratorTemplate: public PseudoGenomeGeneratorBase
    {
        protected:

            // auxiliary structures
            uint_reads_cnt* nextRead = 0;
            uint_read_len* overlap = 0;
            uint_reads_cnt* headRead = 0;
            uint_reads_cnt readsLeft;
            
            bool hasPredecessor(uint_reads_cnt incIdx);
            bool hasSuccessor(uint_reads_cnt incIdx);
            void setReadSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght);

            uint_reads_cnt getHead(uint_reads_cnt idx);
            bool isHeadOf(uint_reads_cnt head, uint_reads_cnt idx);

            uint_pg_len_max pseudoGenomeLength;

            virtual uint_read_len readLength(uint_reads_cnt incIdx) = 0;
            virtual string getReadUpToOverlap(uint_reads_cnt incIdx) = 0;
            virtual uint_reads_cnt readsTotal() = 0;

            virtual ReadsSetProperties* getReadsSetProperties() = 0;
            
            template<typename uint_pg_len, class GeneratedPseudoGenome>
            PseudoGenomeBase* assemblePseudoGenomeFullTemplate();

            template<typename uint_pg_len>
            PseudoGenomeBase* assemblePseudoGenomeTemplate();

            virtual void findOverlappingReads() = 0;
            
            uint_pg_len_max countPseudoGenomeLength();
            uint_reads_cnt countSingles();
            uint_reads_cnt countComponents();
            void quick_stats();

            void init();
            void dispose();
            
        public:

            AbstractOverlapPseudoGenomeGeneratorTemplate() {}
            virtual ~AbstractOverlapPseudoGenomeGeneratorTemplate() {}

            PseudoGenomeBase* generatePseudoGenomeBase();

            bool isPseudoGenomeLengthStandardVirtual();
            bool isPseudoGenomeLengthMaximalVirtual();

    };

}


#endif // ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
