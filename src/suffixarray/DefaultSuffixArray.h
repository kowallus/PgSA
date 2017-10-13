#ifndef DEFAULTSUFFIXARRAY_H_INCLUDED
#define DEFAULTSUFFIXARRAY_H_INCLUDED

#include "SuffixArrayInterface.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "lookuptable/DefaultSuffixArrayLookupTable.h"
#include "../pseudogenome/DefaultPseudoGenome.h"
#include "iterator/OccurrencesIteratorInterface.h"

#include <list>

using namespace PgSAHelpers;

namespace PgSAIndex {

    const string PGSATYPE_DEFAULT = "DEFAULT_PGSA";
    
    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,                
                uchar SA_ELEMENT_SIZE,                // size of sa element in uchar
                uchar POS_START_OFFSET,              // offset of position start offset in uchar
                uint_reads_cnt READSLIST_INDEX_MASK, // LITTLE ENDIAN!
                class ReadsListClass>
    class DefaultSuffixArray: public SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, 
                                DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
                             ,public SuffixArrayBase
                             ,public OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass,
                                DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
    {
        private:

            typedef DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass> ThisDefaultSuffixArrayType;
            typedef DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> PseudoGenome;
            typedef DefaultSuffixArrayLookupTable<uint_read_len, uint_pg_len> SuffixArrayLookupTable;
            typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList;
            typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::ReadsListIteratorClass, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::OccurrencesIteratorClass> OccurrencesIterator;
            
            PseudoGenome* const pseudoGenome;
            ReadsList* readsList; // managed by pseudoGenome
            
            uchar* suffixArray;

            const uint_read_len lookupTableKeyPrefixLength = 11;
            SuffixArrayLookupTable lookupTable;
            
            uchar fixed_min_k;
            
            // HELPER METHODS
            
            const RPGOffset<uint_read_len, uint_reads_cnt> getPositionByAddress(const sa_pos_addr& saPosAddress);

            const char_pg* getSuffixByAddress(const sa_pos_addr& saPosAddress);
            
            const sa_pos_addr saPosIdx2Address(const uint_pg_len posIdx);

            //////////////////////////////////////
            // SUFFIX ARRAY GENERATION ROUTINES //
            //////////////////////////////////////
            
            static const uint_reads_cnt getReadsListIndexByAddress(const sa_pos_addr saPosAddress);

            static const uint_read_len getPosStartOffsetByAddress(const sa_pos_addr saPosAddress);

            static void swapElementsByAddress(const sa_pos_addr saPosAddressFst, const sa_pos_addr saPosAddressSnd);
            static void copyElementsByAddress(const sa_pos_addr saPosAddressDest, const sa_pos_addr saPosAddressSrc);
            
            // auxiliary data for sorting
            static uint_read_len maxReadLength;
            static PseudoGenome* pgStatic;

            static const char_pg* getSuffixStatic(const sa_pos_addr& saPosAddress);

            // TODO: Compare by larger blocks (int or long long int)... needs BIG ENDIAN
            static int pgSuffixesCompare(const void* a, const void* b);

            // The walk around enabling sorting suffix array
            // NOTE: Works only in non-concurent mode (only onfe suffix array at the time)
            
            void prepareUnsortedSA();
            
            void generateSaisPgSA();
            vector<std::ifstream*> saPartSrc;
            int maxPartSize = 0;
            vector<uint_pg_len> currentSaPartPos;
            list<int> saPartsOrder;
            void updateSAPositionQueue(int saGroup);
            
            /////////////////////////////////
            // DUPLICATE FILTER GENERATION //
            /////////////////////////////////
            
            // helper struct (maybe should be local? efficiency tests necessary)
            vector<uint_reads_cnt> readsIdxs;

            // only for constant length reads!
            uint_max markReadsWithDuplicates(uint_max lutIdx);

            void buildReadsWithDuplicatesFilter();
            
            ///////////////////////////////
            // SEARCHING HELPER ROUTINES //
            ///////////////////////////////
            
            typedef int int_kmers_comp;
            
            int_kmers_comp kmerSAComp(const char* kmerPtr, const uint_read_len& kmerLength, const uint_pg_len& suffixIdx, const uint_read_len& lcp);
            
            void kmerRangeBSearch(const char* kmerPtr, const uint_read_len& kmerLength, SARange<uint_pg_len>& range);
            
            uint_max determineSuffixArraySizeInBytesWithGuard(uint_pg_len elementsCount);
            uint_pg_len getSuffixArrayElementsCount(uint_max sizeInBytesWithGuard);
            
        public:
            
            DefaultSuffixArray(PseudoGenome* pseudoGenome, uchar fixed_k = 1);

            DefaultSuffixArray(PseudoGenome* pseudoGenome, std::istream& src);

            void write(ostream& dest);

            ~DefaultSuffixArray();
            
            OccurrencesIterator& getKmerOccurrencesIteratorImpl(const string& kmer);
            
            string getKmerImpl(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength);
            
            void findKmerRangeImpl(const char_pg* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range);
            
            SARange<uint_pg_len> range;
            
            void findOccurrencesOfImpl(const char_pg* kmer, const uint_read_len kmerLength);

            bool moveNextImpl();

            const RPGOffset<uint_read_len, uint_reads_cnt> getPositionImpl(const uint_pg_len posIdx);

            const string getSuffixImpl(const uint_pg_len posIdx, const uint_pg_len length);
            
            ReadsSetInterface<uint_read_len, uint_reads_cnt>* getReadsSetImpl();
            
            string getTypeID();

            const string getDescriptionImpl();
            
    };

    template <  typename uint_read_len,
            typename uint_reads_cnt,
            typename uint_pg_len,
            uchar SA_ELEMENT_SIZE,
            uchar POS_START_OFFSET,
            uint_reads_cnt READSLIST_INDEX_MASK,
            class ReadsListClass>
    uint_read_len DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::maxReadLength;

    template <  typename uint_read_len,
            typename uint_reads_cnt,
            typename uint_pg_len,
            uchar SA_ELEMENT_SIZE,
            uchar POS_START_OFFSET,
            uint_reads_cnt READSLIST_INDEX_MASK,
            class ReadsListClass>
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>*
        DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::pgStatic;
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK>
    using DefaultSuffixArrayOfConstantLengthReadsType = DefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

    template<  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                uchar SA_ELEMENT_SIZE,                // size of sa element in uchar
                uchar POS_START_OFFSET,               // offset of original read index in uchar
                uint_reads_cnt READSLIST_INDEX_MASK // LITTLE INDIAN!
                >
    struct SuffixArrayTraits<DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK>>{
        typedef DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK> SuffixArrayClass; 
        
        typedef typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type ReadsListClass;
        typedef typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass ReadsListIteratorClass;
        typedef SuffixArrayClass OccurrencesIteratorClass;

        typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList; 
        typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, ReadsListIteratorClass, OccurrencesIteratorClass> OccurrencesIterator;        
    };
    
}

#endif // DEFAULTSUFFIXARRAY_H_INCLUDED
