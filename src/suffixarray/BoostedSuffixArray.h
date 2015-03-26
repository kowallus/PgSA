#ifndef BOOSTEDSUFFIXARRAY_H_INCLUDED
#define BOOSTEDSUFFIXARRAY_H_INCLUDED

#include "SuffixArrayInterface.h"
#include "SuffixArrayBase.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "lookuptable/DefaultSuffixArrayLookupTable.h"
#include "../pseudogenome/MultiPackedPseudoGenome.h"
#include "iterator/OccurrencesIteratorInterface.h"

using namespace PgSAHelpers;

namespace PgSAIndex {

    const string PGSATYPE_BOOSTED = "BOOSTED_PGSA";
    
    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,                
                typename uint_pg_element,             // element containing packed symbols in pseudogenome type
                uchar SA_ELEMENT_SIZE,                // size of sa element in uchar
                uchar POS_START_OFFSET,               // offset of position start offset in uchar
                uint_reads_cnt READSLIST_INDEX_MASK,  // LITTLE ENDIAN!
                class ReadsListClass>
    class BoostedSuffixArray: public SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, 
                                BoostedSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
                             ,public SuffixArrayBase
                             ,public OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass,
                                BoostedSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>>
    {
        private:

            typedef BoostedSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, ReadsListClass> ThisDefaultSuffixArrayType;
            typedef MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass> PseudoGenome;
            typedef DefaultSuffixArrayLookupTable<uint_read_len, uint_pg_len> SuffixArrayLookupTable;
            typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList;
            typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::ReadsListIteratorClass, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::OccurrencesIteratorClass> OccurrencesIterator;
                        
            PseudoGenome* const pseudoGenome;
            ReadsList* readsList; // managed by pseudoGenome
            
            uchar* suffixArray;

            uint_read_len lookupTableKeyPrefixLength;
            SuffixArrayLookupTable lookupTable;
            
            // HELPER METHODS
            
            inline static const uint_reads_cnt getReadsListIndexByAddress(const sa_pos_addr saPosAddress) {
                return (*((uint_reads_cnt*) (saPosAddress))) & READSLIST_INDEX_MASK;
            };

            inline static const uint_read_len getPosStartOffsetByAddress(const sa_pos_addr saPosAddress) {
                return *((uint_read_len*) ((uchar*) saPosAddress + POS_START_OFFSET));
            };
            
            inline const RPGOffset<uint_read_len, uint_reads_cnt> getPositionByAddress(const sa_pos_addr& saPosAddress) {
                return { getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress) };
            }

            inline const uint_pg_element* getRawSuffixByAddress(const sa_pos_addr& saPosAddress) {
                return this->pseudoGenome->getRawSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
            };

            inline const string getSuffixByAddress(const sa_pos_addr& saPosAddress, const uint_pg_len length) {
                return this->pseudoGenome->getSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress), length);
            };
            
            inline const sa_pos_addr saPosIdx2Address(const uint_pg_len posIdx) {
                return (sa_pos_addr) (suffixArray + ((uint_max) posIdx) * SA_ELEMENT_SIZE);
            }

            ///////////////////////////////
            // SEARCHING HELPER ROUTINES //
            ///////////////////////////////
            
            typedef int int_kmers_comp;
            
            inline int_kmers_comp kmerSAComp(const uint_pg_element* kmerPtr, const uint_read_len& directCompareElementsCount, const uchar& symbolsLeftToCompare, const uint_pg_len& suffixIdx) {
                const uint_pg_element* suffixPtr = this->getRawSuffix(suffixIdx);

                uint_read_len i = 0;
                while (i++ < directCompareElementsCount) {
                    if (*kmerPtr > *suffixPtr)
                        return 1;
                    if (*kmerPtr++ < *suffixPtr++)
                        return -1;
                }
                if (symbolsLeftToCompare > 0) {
                    uint_pg_element suffixLastValue = this->pseudoGenome->getSymbolsPacker()->clearSuffix(*suffixPtr, symbolsLeftToCompare);
                    if (*kmerPtr > suffixLastValue)
                        return 1;
                    if (*kmerPtr < suffixLastValue)
                        return -1;
                }
                
                return 0;
            }
            
            uint_pg_element packedKmer[USHRT_MAX];
            
            inline void kmerRangeBSearch(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
                this->pseudoGenome->getSymbolsPacker()->packSequence(kmer, kmerLength, packedKmer);
                const uint_read_len directCompareElementsCount = divideBySmallInteger(kmerLength, this->pseudoGenome->getSymbolsPerElement());
                const uchar symbolsLeftToCompare = moduloBySmallInteger(kmerLength, this->pseudoGenome->getSymbolsPerElement(), directCompareElementsCount);  
                
                uint_pg_len mIdx;
                int_kmers_comp cmpRes;

                // lower_bound search
                uint_pg_len lIdx = range.start;
                uint_pg_len rIdx = range.end;

                // (lIdx - 1)           - smaller 
                // (rIdx + 1)           - not smaller
                // range.start - 1      - not greater
                // range.end            - greater 
                
                while (lIdx < rIdx) {
                    mIdx = (lIdx + rIdx) / 2;
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, symbolsLeftToCompare, mIdx);
                    if (cmpRes > 0) {
                        range.start = mIdx + 1;
                        lIdx = mIdx + 1;
                    } else if (cmpRes < 0) {
                        range.end = mIdx;
                        rIdx = mIdx - 1;
                    } else {
                        if(range.start < mIdx)
                            range.start = mIdx + 1;
                        rIdx = mIdx;
                    }
                }

                if (lIdx == rIdx && lIdx == range.start) {
                    // start might be smaller then kmer
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, symbolsLeftToCompare, lIdx);
                    if (cmpRes > 0) 
                        range.start = ++lIdx;
                    else if (cmpRes < 0) {
                        range.end = lIdx;
                        rIdx = lIdx - 1;
                    } else {
                        if(range.start < lIdx)
                            range.start = lIdx + 1;
                        rIdx = lIdx;
                    } 
                }
                
                uint_pg_len start = lIdx;
                
                // upper_bound search
                lIdx = range.start;
                rIdx = range.end;
  
                while (lIdx < rIdx) {
                    mIdx = (lIdx + rIdx) / 2;
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, symbolsLeftToCompare, mIdx);
                    if (cmpRes > 0)
                        lIdx = mIdx + 1;
                    else if (cmpRes < 0)
                        rIdx = mIdx;
                    else
                        lIdx = mIdx + 1;
                }

                range = {start, lIdx};
            }
            
            uint_max getSuffixArraySizeInBytesWithGuard(PseudoGenome* pseudoGenome) {
                return sizeof(uchar) * ((pseudoGenome->getLength()) + 2) * SA_ELEMENT_SIZE;
            }
            
        public:
            
            BoostedSuffixArray(PseudoGenome* pseudoGenome, std::istream& src)
            :   SuffixArrayBase(pseudoGenome->getLength(), getSuffixArraySizeInBytesWithGuard(pseudoGenome), pseudoGenome),
                OccurrencesIterator(
                            *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
                pseudoGenome(pseudoGenome),
                readsList(pseudoGenome->getReadsList()),
                lookupTable(11, pseudoGenome->getReadsSetProperties())
            {
                uint_max arraySize;
                src >> arraySize;
                src.get(); // '/n'

                if (arraySize != this->getSizeInBytes())
                    cout << "WARNING: wrong size of suffixarray.";

                suffixArray = (uchar*) PgSAHelpers::readArray(src, this->getSizeInBytes() * sizeof(uchar));
                src.get(); // '/n'

                this->lookupTable.read(src);

            };

            void write(ostream& dest) {
            };
            
            ~BoostedSuffixArray() {
                delete[] (this->suffixArray);
                delete(pseudoGenome);
            };
            
            inline OccurrencesIterator& getKmerOccurrencesIteratorImpl(const string& kmer) {
                this->findOccurrencesOf(kmer.data(), kmer.length());
                return *this;
            };
            
            inline OccurrencesIterator& getKmerOccurrencesIteratorImpl(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength) {
                const char_pg* kmer = this->pseudoGenome->getSuffixPtrByPosition(originalIdx, pos);
                this->findOccurrencesOf(kmer, kmerLength);
                return *this;
            };
            
            inline void findKmerRangeImpl(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
                if (!lookupTable.findSARange(kmer, kmerLength, range))
                    kmerRangeBSearch(kmer, kmerLength, range);
            }
            
            SARange<uint_pg_len> range;
                        
            inline void findOccurrencesOfImpl(const char* kmer, const uint_read_len& kmerLength) {
                this->findKmerRange(kmer, kmerLength, range);
                this->readsIterator.initIteration(kmerLength);
                if (range.start != range.end)
                    this->readsIterator.setIterationPosition(this->getPosition(range.start), 0);
            }
            
            inline bool moveNextImpl() {
                if (this->readsIterator.moveNext())
                    return true;
                
                while (++(range.start) < range.end) {
                    this->readsIterator.setIterationPosition(this->getPosition(range.start));
                    if (this->readsIterator.moveNext())
                        return true;
                }

                return false;
            };

            inline const RPGOffset<uint_read_len, uint_reads_cnt> getPositionImpl(const uint_pg_len posIdx) {
                return getPositionByAddress(saPosIdx2Address(posIdx));
            };

            inline const uint_pg_element* getRawSuffix(const uint_pg_len posIdx) {
                return getRawSuffixByAddress(saPosIdx2Address(posIdx));
            };
            
            inline const string getSuffixImpl(const uint_pg_len posIdx, const uint_pg_len length) {
                return getSuffixByAddress(saPosIdx2Address(posIdx), length);
            };
            
            ReadsSetInterface<uint_read_len, uint_reads_cnt>* getReadsSetImpl() {
                return this->pseudoGenome;
            }
            
            string getTypeID() { return PGSATYPE_BOOSTED; };

            inline const string getDescriptionImpl() { 
                string desc;
                desc = desc + "Boosted SA (PgSymbolsPerElement = " + toString(this->pseudoGenome->getSymbolsPerElement());
                uint_max size = sizeof(this) + this->getSizeInBytes() + 
                        sizeof(this->lookupTable) + this->lookupTable.getLookupTableLengthWithGuard() * sizeof(uint_pg_len) +
                        sizeof(this->pseudoGenome) + this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement() * this->pseudoGenome->getSymbolsPerElement() +
                        sizeof(this->readsList) + (uint_max) this->readsList->getReadsCount() * this->readsList->getListElementSize();
                desc = desc + ")\tTOTAL (Pg+RL+SA+LT) " + toMB(size, 2) + " MB\n" +
                        + "Pg: " + toMB(this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement() * this->pseudoGenome->getSymbolsPerElement(), 2) + " MB\t" 
                        + "RL: " + toMB((uint_max) this->readsList->getReadsCount() * this->readsList->getListElementSize(), 2) + " MB\t" 
                        + "SA: " + toMB(this->getSizeInBytes(), 2) + " MB\t"
                        + "LT: " + toMB(this->lookupTable.getLookupTableLengthWithGuard() * sizeof(uint_pg_len), 2) + " MB\n";
                desc = desc + toString(this->readsList->getReadsCount()) + " "; 
                if (this->pseudoGenome->isReadLengthConstant())
                    desc = desc + "constant";
                else 
                    desc = desc + "variable";
                desc = desc + " length reads\t max length: " + toString(this->pseudoGenome->maxReadLength()) 
                        + "\t Pg length: " + toString(this->pseudoGenome->getLength()) + "\n";
                return desc;
            };
            
    };
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK>
    using BoostedSuffixArrayOfConstantLengthReadsType = BoostedSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                typename uint_pg_element,             
                uchar SA_ELEMENT_SIZE,                
                uchar POS_START_OFFSET,                    
                uint_reads_cnt READSLIST_INDEX_MASK
                >
    struct SuffixArrayTraits<BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK>>{
        typedef BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, SA_ELEMENT_SIZE, POS_START_OFFSET, READSLIST_INDEX_MASK> SuffixArrayClass; 
        
        typedef typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type ReadsListClass;
        typedef typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass ReadsListIteratorClass;
        typedef SuffixArrayClass OccurrencesIteratorClass;

        typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList; 
        typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, ReadsListIteratorClass, OccurrencesIteratorClass> OccurrencesIterator;        
    };       
    
}

#endif // BOOSTEDSUFFIXARRAY_H_INCLUDED