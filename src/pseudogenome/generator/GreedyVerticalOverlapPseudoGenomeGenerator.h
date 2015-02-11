#ifndef GREEDYVERTICALOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define GREEDYVERTICALOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "../../index/cache/DefaultCountQueriesCache.h"
#include <algorithm>
#include <set>

using namespace PgSAReadsSet;

namespace PgSAIndex {

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead);

    template < typename uint_read_len, typename uint_reads_cnt >
    class GreedyVerticalOverlapGeneratorTemplate: public PseudoGenomeGeneratorBase
    {
        private:

            // auxiliary structures
            uint_reads_cnt* prevRead = 0;
            uint_reads_cnt* nextRead = 0;
            uint_read_len* overlap = 0;
            uint_reads_cnt* headRead = 0;
            uint_reads_cnt readsLeft;
            uint_reads_cnt dupsTotal;

            bool hasPredecessor(uint_reads_cnt incIdx) { return prevRead[incIdx]; };
            bool hasSuccessor(uint_reads_cnt incIdx) { return nextRead[incIdx]; };
            void setReadSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght);

            uint_reads_cnt getHead(uint_reads_cnt idx);
            bool isHeadOf(uint_reads_cnt head, uint_reads_cnt idx) { return getHead(idx) == head; }

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

            uint_pg_len_max pseudoGenomeLength;

            DefaultReadsSet* orgReadsSet = 0;
            uint_read_len readLength(uint_reads_cnt incIdx) { return orgReadsSet->readLength(incIdx - 1); }
            string getRead(uint_reads_cnt incIdx) { return orgReadsSet->getRead(incIdx - 1); }
            uint_reads_cnt readsTotal() { return orgReadsSet->readsCount(); }

            void initOverlapping();
            void populateReadsSet();
            uint_reads_cnt matchInReadsSet(uint_reads_cnt suffixIdx, uint_read_len overlap);
            void overlapMultisetVertical(uint_read_len, uint_read_len);

            template<typename uint_pg_len, class GeneratedPseudoGenome>
            PseudoGenomeBase* assemblePseudoGenomeFullTemplate() {
                clock_checkpoint();
                
                GeneratedPseudoGenome* genPG =
                    new GeneratedPseudoGenome(this->pseudoGenomeLength, this->orgReadsSet->getReadsSetProperties());

                for(uint_reads_cnt i = 1; i <= readsTotal(); i++) {
                    uint_reads_cnt idx = i;
                    if (prevRead[idx] == 0)
                        do {
                            genPG->append(getRead(idx), readLength(idx), overlap[idx], idx - 1);
                            idx = nextRead[idx];
                        } while (idx != 0);
                }

                genPG->validate();
                cout << "Pseudogenome assembled in " << clock_millis() << " msec\n\n";
                
                // TODO: seperate generation of cache
                genPG->setCountQueriesCache(new GeneratedCountQueriesCache(genPG, genPG->getReadsSetProperties()));
                
                return genPG;
            };

            template<typename uint_pg_len>
            PseudoGenomeBase* assemblePseudoGenomeTemplate() {
                if (this->orgReadsSet->isReadLengthConstant())
                    return assemblePseudoGenomeFullTemplate<uint_pg_len, GeneratedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>>();

                cout << "ERROR: Unsuported variable reads length!";
                return 0;
            }

            uint_pg_len_max countPseudoGenomeLength();
            uint_reads_cnt countSingles();
            uint_reads_cnt countComponents();
            void quick_stats();

            void init();
            void dispose();
            
        public:

            GreedyVerticalOverlapGeneratorTemplate(DefaultReadsSet* readsSet);
            virtual ~GreedyVerticalOverlapGeneratorTemplate();

            PseudoGenomeBase* generatePseudoGenomeBase();

            bool isPseudoGenomeLengthStandardVirtual() { return isPGLengthStd(pseudoGenomeLength); };
            bool isPseudoGenomeLengthMaximalVirtual() { return isPGLengthMax(pseudoGenomeLength); };

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
