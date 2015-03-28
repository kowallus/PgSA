#ifndef SPARSESUFFIXARRAY_H_INCLUDED
#define SPARSESUFFIXARRAY_H_INCLUDED

#include "SuffixArrayInterface.h"
#include "SparseSuffixArrayBase.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "lookuptable/DefaultSuffixArrayLookupTable.h"
#include "../pseudogenome/PackedPseudoGenome.h"
#include "iterator/OccurrencesIteratorInterface.h"
#include "../sais/sais.h"

#include <list>

using namespace PgSAHelpers;

namespace PgSAIndex {

    const string PGSATYPE_SPARSE = "SPARSE_PGSA";
    
    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,                
                typename uint_pg_element,             // element containing packed symbols in pseudogenome type
                typename uint_skipped_element,        // element containing packed skipped symbols type
                uchar SA_ELEMENT_SIZE,                // size of sa element in uchar
                uchar POS_START_OFFSET,               // offset of position start offset in uchar
                uchar SKIPPED_OFFSET,        // offset containing packed skipped symbols
                uint_reads_cnt READSLIST_INDEX_MASK,  // LITTLE ENDIAN!
                class ReadsListClass>
    class SparseSuffixArray: public SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, 
                                SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
                             ,public SparseSuffixArrayBase
                             ,public OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass,
                                SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
    {
        private:

            typedef SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass> ThisDefaultSuffixArrayType;
            typedef PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass> PseudoGenome;
            typedef DefaultSuffixArrayLookupTable<uint_read_len, uint_pg_len> SuffixArrayLookupTable;
            typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList;
            typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::ReadsListIteratorClass, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::OccurrencesIteratorClass> OccurrencesIterator;
                        
            PseudoGenome* const pseudoGenome;
            ReadsList* readsList; // managed by pseudoGenome
            
            uchar* suffixArray;
            SymbolsPackingFacility<uint_skipped_element>* sPacker;

            uchar checkLastElementFromShift;
            
            const uint_read_len lookupTableKeyPrefixLength = 11;
            SuffixArrayLookupTable lookupTable;
            
            uchar fixed_min_k;
            
            // HELPER METHODS
            
            const RPGOffset<uint_read_len, uint_reads_cnt> getPositionByAddress(const sa_pos_addr& saPosAddress);

            const uint_pg_element* getRawSuffixByAddress(const sa_pos_addr& saPosAddress);

            const string getSuffixByAddress(const sa_pos_addr& saPosAddress, const uint_pg_len length);
            
            const sa_pos_addr saPosIdx2Address(const uint_pg_len posIdx);
            
            static void swapElementsByAddress(const sa_pos_addr saPosAddressFst, const sa_pos_addr saPosAddressSnd);
            static void copyElementsByAddress(const sa_pos_addr saPosAddressDest, const sa_pos_addr saPosAddressSrc);

            //////////////////////////////////////
            // SUFFIX ARRAY GENERATION ROUTINES //
            //////////////////////////////////////
            
            static const uint_reads_cnt getReadsListIndexByAddress(const sa_pos_addr saPosAddress);

            static const uint_read_len getPosStartOffsetByAddress(const sa_pos_addr saPosAddress);
            
            static const uint_skipped_element getSkippedElementByAddress(const sa_pos_addr saPosAddress);

            // auxiliary data for sorting
            static uint_read_len maxReadLengthRaw;
            static PseudoGenome* pgStatic;

            static const uint_pg_element* getRawSuffixStatic(const sa_pos_addr& saPosAddress);

            // TODO: Compare by larger blocks (int or long long int)... needs BIG ENDIAN
            static int pgSuffixesCompare(const void* a, const void* b);

            uint_pg_len elementsCountWithoutGuard = 0;
            
            void prepareUnsortedSA();
            
            void generatePgSA();

            vector<std::ifstream*> saPartSrc;
            int maxPartSize = 0;
            vector<uint_pg_len> currentSaPartPos;
            list<int> saPartsOrder;
            
            void updateSAPositionQueue(int saGroup);
            
            void generateSaisPgSA();

            void buildReadsWithDuplicatesFilter();
            
            uint_max markReadsWithDuplicates(string kmer);
            
            ///////////////////////////////
            // SEARCHING HELPER ROUTINES //
            ///////////////////////////////
            
            typedef int int_kmers_comp;
            
            int_kmers_comp kmerSAComp(const uint_pg_element* kmerPtr, const uint_read_len& directCompareElementsCount, const uint_read_len& omitElementsCount, const uchar& symbolsLeftToCompare, const uint_pg_len& suffixIdx);
            
            uint_pg_element packedKmer[USHRT_MAX];
            
            void kmerRangeBSearch(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range);

            const uint_max determineElementsCount(PseudoGenome* pseudoGenome);
            
            uint_max getSuffixArraySizeInBytesWithGuard(uint_pg_len elementsCount);
            uint_pg_len getSuffixArrayElementsCount(uint_max sizeInBytesWithGuard);
            
            void determineCheckLastElementFromShift(uint_read_len& kmerLength);
            
        public:
            
            SparseSuffixArray(PseudoGenome* pseudoGenome, int fixed_min_k = 1);

            SparseSuffixArray(PseudoGenome* pseudoGenome, std::istream& src);

            void write(ostream& dest);
            
            ~SparseSuffixArray();
            
            OccurrencesIterator& getKmerOccurrencesIteratorImpl(const string& kmer);
            
            char_pg kmerByPos[USHRT_MAX];
            
            OccurrencesIterator& getKmerOccurrencesIteratorImpl(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength);
            
            void findKmerRangeImpl(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range);
            
            const char* kmer;
            uint_read_len kmerLength;
            SARange<uint_pg_len> range;
            uchar shift;
            
            bool matchPrefix(const sa_pos_addr& saPosAddress);
            
            void findOccurrencesOfImpl(const char* kmer, const uint_read_len& kmerLength);
            
            void findOccurrences(uchar shift);

            bool moveNextImpl();

            const RPGOffset<uint_read_len, uint_reads_cnt> getPositionImpl(const uint_pg_len posIdx);

            const uint_pg_element* getRawSuffix(const uint_pg_len posIdx);
            
            const string getSuffixImpl(const uint_pg_len posIdx, const uint_pg_len length);
            
            ReadsSetInterface<uint_read_len, uint_reads_cnt>* getReadsSetImpl();
            
            string getTypeID();

            const string getDescriptionImpl();
            
    };

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                typename uint_pg_element,             
                typename uint_skipped_element,  
                uchar SA_ELEMENT_SIZE,       
                uchar POS_START_OFFSET,               
                uchar SKIPPED_OFFSET,        
                uint_reads_cnt READSLIST_INDEX_MASK,        
                class ReadsListClass>
    uint_read_len SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::maxReadLengthRaw;

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                typename uint_pg_element,        
                typename uint_skipped_element,
                uchar SA_ELEMENT_SIZE,
                uchar POS_START_OFFSET,               
                uchar SKIPPED_OFFSET,        
                uint_reads_cnt READSLIST_INDEX_MASK,  
                class ReadsListClass>
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::pgStatic;
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element , uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK>
    using SparseSuffixArrayOfConstantLengthReadsType = SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                typename uint_pg_element,             
                typename uint_skipped_element,
                uchar SA_ELEMENT_SIZE,                
                uchar POS_START_OFFSET,               
                uchar SKIPPED_OFFSET,        
                uint_reads_cnt READSLIST_INDEX_MASK
                >
    struct SuffixArrayTraits<SparseSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK>>{
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK> SuffixArrayClass; 
        
        typedef typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type ReadsListClass;
        typedef typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass ReadsListIteratorClass;
        typedef SuffixArrayClass OccurrencesIteratorClass;

        typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList; 
        typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, ReadsListIteratorClass, OccurrencesIteratorClass> OccurrencesIterator;        
    };       
    
}

#endif // SPARSESUFFIXARRAY_H_INCLUDED