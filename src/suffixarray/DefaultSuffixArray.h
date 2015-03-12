#ifndef DEFAULTSUFFIXARRAY_H_INCLUDED
#define DEFAULTSUFFIXARRAY_H_INCLUDED

#include "SuffixArrayInterface.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "lookuptable/DefaultSuffixArrayLookupTable.h"
#include "../pseudogenome/DefaultPseudoGenome.h"
#include "iterator/OccurrencesIteratorInterface.h"

#include <byteswap.h>

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
            typedef DefaultSuffixArrayLookupTable<uint_read_len, uint_pg_len, ThisDefaultSuffixArrayType> SuffixArrayLookupTable;
            typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList;
            typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::ReadsListIteratorClass, typename SuffixArrayTraits<ThisDefaultSuffixArrayType>::OccurrencesIteratorClass> OccurrencesIterator;
            
            PseudoGenome* const pseudoGenome;
            ReadsList* readsList; // managed by pseudoGenome
            
            uchar* suffixArray;

            const uint_read_len lookupTableKeyPrefixLength = 11;
            SuffixArrayLookupTable lookupTable;
            
            // HELPER METHODS
            
            inline const RPGOffset<uint_read_len, uint_reads_cnt> getPositionByAddress(const sa_pos_addr& saPosAddress) {
                return { getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress) };
            }

            inline const char_pg* getSuffixByAddress(const sa_pos_addr& saPosAddress) {
                return this->pseudoGenome->getSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
            };
            
            inline const sa_pos_addr saPosIdx2Address(const uint_pg_len posIdx) {
                return (sa_pos_addr) (suffixArray + ((uint_max) posIdx) * SA_ELEMENT_SIZE);
            }

            //////////////////////////////////////
            // SUFFIX ARRAY GENERATION ROUTINES //
            //////////////////////////////////////
            
            inline static const uint_reads_cnt getReadsListIndexByAddress(const sa_pos_addr saPosAddress) {
                return (*((uint_reads_cnt*) (saPosAddress))) & READSLIST_INDEX_MASK;
            };

            inline static const uint_read_len getPosStartOffsetByAddress(const sa_pos_addr saPosAddress) {
                return *((uint_read_len*) ((uchar*) saPosAddress + POS_START_OFFSET));
            };

            // auxiliary data for sorting
            static uint_read_len maxReadLength;
            static PseudoGenome* pgStatic;

            inline static const char_pg* getSuffixStatic(const sa_pos_addr& saPosAddress) {
                return pgStatic->getSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
            }

            // TODO: Compare by larger blocks (int or long long int)... needs BIG ENDIAN
            static int pgSuffixesCompare(const void* a, const void* b) {
                const char_pg* readA = getSuffixStatic((sa_pos_addr) a);
                const char_pg* readB = getSuffixStatic((sa_pos_addr) b);

                int i = 0;
                while (i++ < maxReadLength) {
                    if (*readA > *readB)
                        return 1;
                    if (*readA++ < *readB++)
                        return -1;
                }
                return 0;
            }

            // The walk around enabling sorting suffix array
            // NOTE: Works only in non-concurent mode (only onfe suffix array at the time)
            void generatePgSA() {
                clock_checkpoint();
                
                const uchar* curSAPos = suffixArray;

                uint_reads_cnt readsListIndex = 0;

                for(uint_pg_len pgPos = 0; pgPos < this->elementsCount; pgPos++) {
                    while (readsList->getReadPosition(readsListIndex + 1) <= pgPos)
                        readsListIndex++ ;
                    *((uint_reads_cnt*) curSAPos) = readsListIndex;
                    *((uint_read_len*) (curSAPos + POS_START_OFFSET)) = pgPos - readsList->getReadPosition(readsListIndex);
                    curSAPos += SA_ELEMENT_SIZE;
                }

                if (curSAPos != suffixArray + this->pseudoGenome->getLength() * (uint_max) SA_ELEMENT_SIZE )
                    cout << "WARNING: SA generation failed: " << (int) (curSAPos - suffixArray) / SA_ELEMENT_SIZE << " elements instead of " << this->pseudoGenome->getLength() << "\n";

                pgStatic = this->pseudoGenome;
                maxReadLength = this->pseudoGenome->maxReadLength();

                qsort(suffixArray, this->elementsCount, sizeof(uchar) * SA_ELEMENT_SIZE, this->pgSuffixesCompare);                
                
                cout << "SA generation time " << clock_millis() << " msec!\n";
            }

            /////////////////////////////////
            // DUPLICATE FILTER GENERATION //
            /////////////////////////////////
            
            // helper struct (maybe should be local? efficiency tests necessary)
            vector<uint_reads_cnt> readsIdxs;

            // only for constant length reads!
            uint_max markReadsWithDuplicates(uint_max lutIdx) {

                readsIdxs.clear();
                uint_pg_len start = lookupTable.getRawValue(lutIdx);
                uint_pg_len stop =  lookupTable.getRawValue(lutIdx + 1);
                const uint_read_len guardOffset = readsList->getMaxReadLength() - readsList->getDuplicateFilterKmerLength() ;

                for (uint_pg_len i = start; i < stop; i++) {
                    const sa_pos_addr saPosAddress = this->saPosIdx2Address(i);
                    uint_reads_cnt j = this->getReadsListIndexByAddress(saPosAddress);
                    uint_pg_len guard = this->getPosStartOffsetByAddress(saPosAddress) + readsList->getReadPosition(j) - guardOffset;
                    while (readsList->getReadPosition(j) >= guard) {
                        readsList->setOccurOnceFlag(j);
                        readsIdxs.push_back(j);
                        if (j == 0) break;
                        j--;
                    }
                }

                uint_reads_cnt j = 0;
                uint_reads_cnt readsCount = readsIdxs.size();
                for(uint_reads_cnt i = 0; i < readsCount; i++) {
                    if (!readsList->hasOccurOnceFlag(readsIdxs[i]) && !readsList->hasDuplicateFilterFlag(readsIdxs[i])) {
                        readsList->setDuplicateFilterFlag(readsIdxs[i]);
                        j++;
                    }
                    else
                        readsList->clearOccurFlags(readsIdxs[i]);
                }

//                if (j > 0)
//                    cout << "Found " << j << " reads with duplicates for lutIdx " << lutIdx << "\n";

                return j;
            }

            void buildReadsWithDuplicatesFilter() {
                clock_checkpoint();
                
                readsList->setDuplicateFilterKmerLength(this->lookupTable.getKeyPrefixLength());
                if (readsList->getDuplicateFilterKmerLength() != this->lookupTable.getKeyPrefixLength()) {
                    cout << "Unsupported duplicate filter size " << readsList->getDuplicateFilterKmerLength() << " expected " << this->lookupTable.getKeyPrefixLength() << "!\n";
                    exit(-1);
                }

                uint_max filterCount = 0;
                uint_max lutSize = lookupTable.getLookupTableLengthWithGuard();

                for(uint_max j = 0; j < lutSize - 1; j++)
                    filterCount += markReadsWithDuplicates(j);

                cout << "Found " << filterCount << " reads containing duplicate " << (int) readsList->getDuplicateFilterKmerLength() 
                        << "-mers in " << clock_millis() << " msec!\n";
            }

            ///////////////////////////////
            // SEARCHING HELPER ROUTINES //
            ///////////////////////////////
            
            typedef int int_kmers_comp;
            
            inline int_kmers_comp kmerSAComp(const char* kmerPtr, const uint_read_len& kmerLength, const uint_pg_len& suffixIdx, const uint_read_len& lcp) {
                kmerPtr += lcp;
                const char* suffixPtr = getSuffixByAddress(saPosIdx2Address(suffixIdx)) + lcp;

                uint_read_len i = lcp;

                while (kmerLength - i >= 2) {
                    int cmp = bswap_16(*(uint16_t*)kmerPtr) - bswap_16(*(uint16_t*)suffixPtr);
                    if (cmp != 0)
                        break;
                    kmerPtr += 2;
                    suffixPtr += 2;
                    i += 2;
                }

                while (i < kmerLength) {
                    i++;
                    int cmp = *kmerPtr++ - *suffixPtr++;
                    if (cmp > 0)
                        return i;
                    if (cmp < 0)
                        return -i;
                }

                return 0;
                
            }
            
            inline void kmerRangeBSearch(const char* kmerPtr, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
                
                uint_pg_len mIdx;
                int_kmers_comp cmpRes;

                uint_read_len lcp = lookupTableKeyPrefixLength;
                
                uint_read_len lcp_l = lcp;
                uint_read_len lcp_r = lcp;
                uint_read_len lcp_beg = lcp;
                uint_read_len lcp_end = lcp;
                
                // lower_bound search
                uint_pg_len lIdx = range.start;
                uint_pg_len rIdx = range.end;

                while (lIdx < rIdx) {
                    mIdx = (lIdx + rIdx) / 2;
                    cmpRes = kmerSAComp(kmerPtr, kmerLength, mIdx, lcp);
              
                    if (cmpRes > 0) {
                        range.start = mIdx + 1;
                        lIdx = mIdx + 1;
                        
                        if ((lcp_l = cmpRes - 1) <= lcp_r)
                            lcp = lcp_l;
                        else if (lcp_r > lcp) 
                            lcp = lcp_r;
                        lcp_beg = lcp_l;
                    } else if (cmpRes < 0) {
                        range.end = mIdx;
                        rIdx = mIdx - 1;
                        
                        if ((lcp_r = -cmpRes - 1) <= lcp_l)                            
                            lcp = lcp_r;
                        else if (lcp_l > lcp)
                            lcp = lcp_l;
                        lcp_end = lcp_r;
                    } else {
                        if(range.start < mIdx)
                            range.start = mIdx + 1;
                        rIdx = mIdx;
                        
                        lcp_r = kmerLength;
                        lcp_beg = lcp_r;
                    }
                }

                if (lIdx == rIdx && lIdx == range.start) {
                    // start might be smaller then kmer
                    cmpRes = kmerSAComp(kmerPtr, kmerLength, lIdx, lcp);
                    
                    if (cmpRes > 0) {
                        range.start = ++lIdx;
                        lcp_beg = lcp_l;
                    } else if (cmpRes < 0) {
                        range.end = lIdx;
                        lcp_end = lcp_l;
                    } else {
                        if(range.start < lIdx)
                            range.start = lIdx + 1;
                        lcp_beg = lcp_l;
                    } 
                }
                
                uint_pg_len start = lIdx;

                // upper_bound search
                lIdx = range.start;
                rIdx = range.end;
                lcp_l = lcp_beg;
                lcp_r = lcp_end;
                lcp = lcp_l < lcp_r?lcp_l:lcp_r;
                
                while (lIdx < rIdx) {
                    mIdx = (lIdx + rIdx) / 2;
                    cmpRes = kmerSAComp(kmerPtr, kmerLength, mIdx, lcp);
                    
                    if (cmpRes > 0) {
                        lIdx = mIdx + 1;
                      
                        if ((lcp_l = cmpRes - 1) <= lcp_r)
                            lcp = lcp_l;
                        else if (lcp_r > lcp) 
                            lcp = lcp_r;
                    } else if (cmpRes < 0) {
                        if ((lcp_r = -cmpRes - 1) <= lcp_l)
                            lcp = lcp_r;
                        else if (lcp_l > lcp) 
                            lcp = lcp_l;
                    } else {
                        lIdx = mIdx + 1;
                        
                        lcp_l = kmerLength;
                    }
                }

                range = {start, lIdx};
            }
            
            uint_max getSuffixArraySizeInBytesWithGuard(PseudoGenome* pseudoGenome) {
                return sizeof(uchar) * ((pseudoGenome->getLength()) + 2) * SA_ELEMENT_SIZE;
            }
            
        public:
            
            DefaultSuffixArray(PseudoGenome* pseudoGenome)
            :   // up to 2 additional bytes to avoid overflowing during casting to uint_reads_cnt
                SuffixArrayBase(pseudoGenome->getLength(), getSuffixArraySizeInBytesWithGuard(pseudoGenome), pseudoGenome),
                OccurrencesIterator(
                            *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
                pseudoGenome(pseudoGenome),
                readsList(pseudoGenome->getReadsList()),
                suffixArray(new uchar[this->getSizeInBytes()]),
                lookupTable(lookupTableKeyPrefixLength, pseudoGenome->getReadsSetProperties())
            {
                //TODO: check if setting 0 is necessary.... 
                suffixArray[this->getSizeInBytes() - 2] = 0;
                suffixArray[this->getSizeInBytes() - 1] = 0;
                
 //               cout << "memcheck 3.....\n";
//                cin.ignore(1);
                generatePgSA();
                this->lookupTable.generate(this, this->getElementsCount());
//                cout << "memcheck 4.....\n";
//                cin.ignore(1);
                buildReadsWithDuplicatesFilter();
//                cout << "memcheck 5.....\n";
//                cin.ignore(1);
            };

            DefaultSuffixArray(PseudoGenome* pseudoGenome, std::istream& src)
            :   SuffixArrayBase(pseudoGenome->getLength(), getSuffixArraySizeInBytesWithGuard(pseudoGenome), pseudoGenome),
                OccurrencesIterator(
                            *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
                pseudoGenome(pseudoGenome),
                readsList(pseudoGenome->getReadsList()),
                lookupTable(lookupTableKeyPrefixLength, pseudoGenome->getReadsSetProperties())
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
                dest << this->getSizeInBytes() << "\n";
                PgSAHelpers::writeArray(dest, suffixArray, this->getSizeInBytes());
                dest << "\n";
                lookupTable.write(dest);
            };

            ~DefaultSuffixArray() {
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
            
            inline void findKmerRangeImpl(const char_pg* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
                if (!lookupTable.findSARange(kmer, kmerLength, range))
                    kmerRangeBSearch(kmer, kmerLength, range);
            }
            
            SARange<uint_pg_len> range;
            
            inline void findOccurrencesOfImpl(const char_pg* kmer, const uint_read_len kmerLength) { 
                this->findKmerRange(kmer, kmerLength, range);
                this->readsIterator.initIteration(kmerLength);              

                if (range.start != range.end)
                    this->readsIterator.setIterationPosition(this->getPosition(range.start));
            };

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

            inline const string getSuffixImpl(const uint_pg_len posIdx, const uint_pg_len length) {
                return string(getSuffixByAddress(saPosIdx2Address(posIdx)), length);
            };
            
            ReadsSetInterface<uint_read_len, uint_reads_cnt>* getReadsSetImpl() {
                return this->pseudoGenome;
            }
            
            string getTypeID() { return PGSATYPE_DEFAULT; };

            inline const string getDescriptionImpl() { 
                string desc;
                uint_max size = sizeof(this) + this->getSizeInBytes() + 
                        sizeof(this->lookupTable) + this->lookupTable.getLookupTableLengthWithGuard() * sizeof(uint_pg_len) +
                        sizeof(this->pseudoGenome) + this->pseudoGenome->getLengthWithGuard() * sizeof(char_pg) +
                        sizeof(this->readsList) + (uint_max) this->readsList->getReadsCount() * this->readsList->getListElementSize();
                desc = desc + "Standard SA\t TOTAL (Pg+RL+SA+LT) " + toMB(size, 2) + " MB\n"
                        + "Pg: " + toMB(this->pseudoGenome->getLengthWithGuard() * sizeof(char_pg), 2) + " MB\t" 
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
