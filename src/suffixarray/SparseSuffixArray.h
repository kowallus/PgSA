#ifndef SPARSESUFFIXARRAY_H_INCLUDED
#define SPARSESUFFIXARRAY_H_INCLUDED

#include "SuffixArrayInterface.h"
#include "SparseSuffixArrayBase.h"
#include "../pseudogenome/PseudoGenomeInterface.h"
#include "lookuptable/DefaultSuffixArrayLookupTable.h"
#include "../pseudogenome/PackedPseudoGenome.h"
#include "iterator/OccurrencesIteratorInterface.h"
#include "../sais/sais.h"

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
            
            // HELPER METHODS
            
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
            
            inline static void swapElementsByAddress(const sa_pos_addr saPosAddressFst, const sa_pos_addr saPosAddressSnd) {
                uchar tmp;
                for(int i = 0; i < SA_ELEMENT_SIZE; i++) {
                    tmp = *(((uchar*) saPosAddressFst)+i);
                    *(((uchar*) saPosAddressFst)+i) = *(((uchar*) saPosAddressSnd)+i);
                    *(((uchar*) saPosAddressSnd)+i) = tmp;
        }
            };
            

            //////////////////////////////////////
            // SUFFIX ARRAY GENERATION ROUTINES //
            //////////////////////////////////////
            
            inline static const uint_reads_cnt getReadsListIndexByAddress(const sa_pos_addr saPosAddress) {
                return (*((uint_reads_cnt*) (saPosAddress))) & READSLIST_INDEX_MASK;
            };

            inline static const uint_read_len getPosStartOffsetByAddress(const sa_pos_addr saPosAddress) {
                return *((uint_read_len*) ((uchar*) saPosAddress + POS_START_OFFSET));
            };
            
            inline static const uint_skipped_element getSkippedElementByAddress(const sa_pos_addr saPosAddress) {
                return *((uint_skipped_element*) ((uchar*) saPosAddress + SKIPPED_OFFSET));
            };

            // auxiliary data for sorting
            static uint_read_len maxReadLengthRaw;
            static PseudoGenome* pgStatic;

            inline static const uint_pg_element* getRawSuffixStatic(const sa_pos_addr& saPosAddress) {
                return pgStatic->getRawSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
            }

            // TODO: Compare by larger blocks (int or long long int)... needs BIG ENDIAN
            static int pgSuffixesCompare(const void* a, const void* b) {
                const uint_pg_element* readA = getRawSuffixStatic((sa_pos_addr) a);
                const uint_pg_element* readB = getRawSuffixStatic((sa_pos_addr) b);

                int i = 0;
                while (i++ < maxReadLengthRaw) {
                    if (*readA > *readB)
                        return 1;
                    if (*readA++ < *readB++)
                        return -1;
                }
                return 0;
            }

            uint_pg_len elementsCountWithoutGuard = 0;
            
            void prepareUnsortedSA() {
                
                const uchar* curSAPos = suffixArray;

                uint_reads_cnt readsListIndex = 0;

                for(uint_pg_len pgPos = 0; pgPos <= pseudoGenome->getLength(); pgPos += pseudoGenome->getSymbolsPerElement()) {
                    uint_pg_len pgSuffixPos = pgPos + skippedSymbolsCount;
                    while ((readsListIndex < readsList->getReadsCount()) && 
                            (readsList->getReadPosition(readsListIndex + 1) <= pgSuffixPos))
                        readsListIndex++ ;
                    *((uint_reads_cnt*) curSAPos) = readsListIndex;
                    *((uint_read_len*) (curSAPos + POS_START_OFFSET)) = pgSuffixPos - readsList->getReadPosition(readsListIndex);
                    *((uint_skipped_element*) (curSAPos + SKIPPED_OFFSET)) = 
                            sPacker->packSymbols(pseudoGenome->getSuffix(pgPos, skippedSymbolsCount).data());
                            
                    curSAPos += SA_ELEMENT_SIZE;
                }

                elementsCountWithoutGuard = this->getPosition(this->elementsCount - 1).readListIndex == this->pseudoGenome->readsCount()?this->elementsCount - 1:this->elementsCount;
                
                if (curSAPos != suffixArray + this->elementsCount * (uint_max) SA_ELEMENT_SIZE )
                    cout << "WARNING: SA generation failed: " << (int) (curSAPos - suffixArray) / SA_ELEMENT_SIZE << " elements instead of " << this->elementsCount << "\n";

            }
            
            void generatePgSA() {
                clock_checkpoint();
                
                prepareUnsortedSA();
                
                qsort(suffixArray, elementsCountWithoutGuard, sizeof(uchar) * SA_ELEMENT_SIZE, this->pgSuffixesCompare);
                
                cout << "SA generation time " << clock_millis() << " msec!\n";
            }

            vector<std::ifstream*> saPartSrc;
            int maxPartSize = 0;
            vector<uint_pg_len> currentSaPartPos;
            list<int> saPartsOrder;
            
            void updateSAPositionQueue(int saGroup) {
                int relativePosition = 0;

                while (saPartSrc[saGroup]->read((char*) &relativePosition, sizeof(int))) {

                    if (relativePosition >= maxPartSize)
                        continue;
                    currentSaPartPos[saGroup] = (uint_pg_len) saGroup * (uint_pg_len) maxPartSize + ((uint_pg_len) relativePosition);
                    if (currentSaPartPos[saGroup] >= elementsCountWithoutGuard)
                        continue;

                    list<int>::reverse_iterator it = saPartsOrder.rbegin();
                    while (true) {
                        if (it == saPartsOrder.rend() || (strcmplcp((const char*) pseudoGenome->getRawSuffix(1 + currentSaPartPos[saGroup]), (const char*) pseudoGenome->getRawSuffix(1 + currentSaPartPos[*it]), maxReadLengthRaw) >= 0)) {
                            saPartsOrder.insert(it.base(), saGroup);
                            return;
                        }
                        it++;
                    }
                }
            }
            
            void generateSaisPgSA() {
                clock_checkpoint();
        
                int maxPartBruttoSize = INT_MAX / (sizeof(int));
                uint_pg_len pgElementsCountWithGuard = pseudoGenome->getElementsCountWithGuard();

                if ((uint_max) maxPartBruttoSize * 6 > pgElementsCountWithGuard * (uint_max) SA_ELEMENT_SIZE)
                    maxPartBruttoSize = pgElementsCountWithGuard * (uint_max) SA_ELEMENT_SIZE / 6;

                maxPartSize = maxPartBruttoSize - maxReadLengthRaw;
                int noOfParts = 0; 

                int* saisSA = (int*) malloc((size_t)(maxPartBruttoSize * sizeof(int)));

                while (noOfParts * (uint_pg_len) maxPartSize < pgElementsCountWithGuard ) {

                    int partSize = pgElementsCountWithGuard - noOfParts * maxPartSize;
                    if (partSize > maxPartBruttoSize)
                        partSize = maxPartBruttoSize;

                    if(sais((const unsigned char*) this->pseudoGenome->getRawSuffix(1 + noOfParts * maxPartSize), saisSA, (int) partSize) != 0) {
                        fprintf(stderr, "Cannot allocate memory.\n");
                        exit(EXIT_FAILURE);
                    }

                    std::ofstream dest("saisSA" + toString(noOfParts++) + ".tmp", std::ios::out | std::ios::binary);
                    PgSAHelpers::writeArray(dest, saisSA, (size_t)(partSize) * sizeof(int));
                    dest.close();
                }

                free((void*)saisSA);

                cout << "SAIS generation time " << clock_millis() << " msec!\n";

                clock_checkpoint();

                prepareUnsortedSA();

                readsList->buildLUT();

                currentSaPartPos.resize(noOfParts);
                for (int i = 0; i < noOfParts; i++) {
                    saPartSrc.push_back(new std::ifstream("saisSA" + toString(i) + ".tmp", std::ifstream::binary));
                    updateSAPositionQueue(i);
                }

                const uchar* curSAPos = suffixArray;
                uint_pg_len saPos;
                for (uint_pg_len i = 0; i < elementsCountWithoutGuard; i++) {
                    int j = saPartsOrder.front();
                    saPos = currentSaPartPos[j];

                    uint_pg_len pgPos = saPos * pseudoGenome->getSymbolsPerElement() + skippedSymbolsCount;
                    
                    if ((uint_pg_len) saPos < i) {
                        uint_reads_cnt readsListIndex = readsList->findFurthestReadContaining(pgPos);
                        *((uint_reads_cnt*) curSAPos) = readsListIndex;
                        *((uint_read_len*) (curSAPos + POS_START_OFFSET)) = pgPos - readsList->getReadPosition(readsListIndex);
                        *((uint_skipped_element*) (curSAPos + SKIPPED_OFFSET)) = 
                            sPacker->packSymbols(pseudoGenome->getSuffix(pgPos - skippedSymbolsCount, skippedSymbolsCount).data());
                    } else
                        swapElementsByAddress((sa_pos_addr*) curSAPos, saPosIdx2Address(saPos));

                    curSAPos += SA_ELEMENT_SIZE;

                    saPartsOrder.pop_front();
                    updateSAPositionQueue(j);
                }

                for (int i = 0; i < noOfParts; i++) {
                    saPartSrc[i]->close();
                    delete(saPartSrc[i]);
                    remove(("saisSA" + toString(noOfParts) + ".tmp").c_str());
                }

                saPartSrc.clear();

                cout << "SA generation time " << clock_millis() << " msec!\n";
            }

            void buildReadsWithDuplicatesFilter() {
                clock_checkpoint();

                readsList->setDuplicateFilterKmerLength(this->lookupTable.getKeyPrefixLength());
                if (readsList->getDuplicateFilterKmerLength() != this->lookupTable.getKeyPrefixLength()) {
                    cout << "Unsupported duplicate filter size " << readsList->getDuplicateFilterKmerLength() << " expected " << this->lookupTable.getKeyPrefixLength() << "!\n";
                    exit(-1);
                }
                uint_max filterCount = 0;
                
                ReadsSetProperties* properties = this->pseudoGenome->getReadsSetProperties();
                
                string kmer(this->lookupTable.getKeyPrefixLength(), properties->symbolsList[0]);
                               
                while (true) {
                    filterCount += markReadsWithDuplicates(kmer);
                    
                    uchar i = this->lookupTable.getKeyPrefixLength();
                    uchar order;
                    do {
                        order = properties->symbolOrder[(unsigned char) kmer[--i]] + 1;
                        if (order == properties->symbolsCount)
                            order = 0;
                        kmer[i] = properties->symbolsList[order];
                    } while (order == 0 && i > 0);
                    if (i == 0 && order == 0)
                        break;
                }
                
                cout << "Found " << filterCount << " reads containing duplicate " << (int) readsList->getDuplicateFilterKmerLength()
                        << "-mers in " << clock_millis() << " msec!\n";
            }
            
            uint_max markReadsWithDuplicates(string kmer) {

                vector<uint_reads_cnt> readsIdxs;
                
                const uint_read_len guardOffset = readsList->getMaxReadLength() - readsList->getDuplicateFilterKmerLength();
                OccurrencesIterator& oit = this->getKmerOccurrencesIterator(kmer);
                
                while(oit.moveNext()) {
                    uint_reads_cnt j = oit.getReadIndex();
                     
                    if (readsList->hasDuplicateFilterFlag(j))
                        continue;
                    if (guardOffset >= oit.getOccurrenceOffset()) {
                        if (!readsList->hasOccurFlag(j))
                            readsIdxs.push_back(j);
                        readsList->setOccurOnceFlag(j);
                    }
                }

                uint_reads_cnt j = 0;
                uint_reads_cnt readsCount = readsIdxs.size();
                for (uint_reads_cnt i = 0; i < readsCount; i++) {
                    if (!readsList->hasOccurOnceFlag(readsIdxs[i]) && !readsList->hasDuplicateFilterFlag(readsIdxs[i])) {
                        readsList->setDuplicateFilterFlag(readsIdxs[i]);
                        j++;               
                    } 
                    readsList->clearOccurFlags(readsIdxs[i]);
                }

                return j;
            }
            
            ///////////////////////////////
            // SEARCHING HELPER ROUTINES //
            ///////////////////////////////
            
            typedef int int_kmers_comp;
            
            inline int_kmers_comp kmerSAComp(const uint_pg_element* kmerPtr, const uint_read_len& directCompareElementsCount, const uint_read_len& omitElementsCount, const uchar& symbolsLeftToCompare, const uint_pg_len& suffixIdx) {
                kmerPtr += omitElementsCount;
                const uint_pg_element* suffixPtr = this->getRawSuffix(suffixIdx) + omitElementsCount;

                uint_read_len i = omitElementsCount;
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
                const uint_read_len omitElementsCount = divideBySmallInteger(lookupTableKeyPrefixLength, this->pseudoGenome->getSymbolsPerElement());
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
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, omitElementsCount, symbolsLeftToCompare, mIdx);
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
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, omitElementsCount, symbolsLeftToCompare, lIdx);
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
                    cmpRes = kmerSAComp(packedKmer, directCompareElementsCount, omitElementsCount, symbolsLeftToCompare, mIdx);
                    if (cmpRes > 0)
                        lIdx = mIdx + 1;
                    else if (cmpRes < 0)
                        rIdx = mIdx;
                    else
                        lIdx = mIdx + 1;
                }

                range = {start, lIdx};
            }

            const uint_max determineElementsCount(PseudoGenome* pseudoGenome) {
                return pseudoGenome->getLength() / pseudoGenome->getSymbolsPerElement() + 1;
            }
            
            const uint_max getSizeInBytesWithGuard(PseudoGenome* pseudoGenome) {
                return sizeof(uchar) * (pseudoGenome->getLength() / pseudoGenome->getSymbolsPerElement() + 2) * SA_ELEMENT_SIZE;
            }
            
            void determineCheckLastElementFromShift(uint_read_len& kmerLength) {
                RPGOffset<uint_read_len, uint_reads_cnt> lastElement = this->getPosition(this->getElementsCount() - 1);
                if (lastElement.readListIndex == pseudoGenome->readsCount())
                    checkLastElementFromShift = lastElement.offset+kmerLength;
                else
                    checkLastElementFromShift = 0;
            }
            
        public:
            
            SparseSuffixArray(PseudoGenome* pseudoGenome)
            :   // up to 2 additional bytes to avoid overflowing during casting to uint_reads_cnt
                SparseSuffixArrayBase(determineElementsCount(pseudoGenome), getSizeInBytesWithGuard(pseudoGenome), pseudoGenome, pseudoGenome->getSymbolsPerElement()-1, sizeof(uint_skipped_element)),
                OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass, ThisDefaultSuffixArrayType>(
                            *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
                pseudoGenome(pseudoGenome),
                readsList(pseudoGenome->getReadsList()),
                suffixArray(new uchar[this->getSizeInBytes()]),
                sPacker(new SymbolsPackingFacility<uint_skipped_element>(pseudoGenome->getReadsSetProperties(), pseudoGenome->getSymbolsPerElement()-1)),
                lookupTable(lookupTableKeyPrefixLength, pseudoGenome->getReadsSetProperties())
            {
                //TODO: check if setting 0 is necessary.... 
                suffixArray[this->getSizeInBytes() - 2] = 0;
                suffixArray[this->getSizeInBytes() - 1] = 0;
                
                pgStatic = this->pseudoGenome;
                maxReadLengthRaw = (this->pseudoGenome->maxReadLength() + this->pseudoGenome->getSymbolsPerElement() - 1) / this->pseudoGenome->getSymbolsPerElement();
                
                this->lookupTable.generateFromPg(this->pseudoGenome, pseudoGenome->getSymbolsPerElement(), skippedSymbolsCount);
                
                if ((sizeof(uint_pg_element) == sizeof(char)))
                    generateSaisPgSA();
                else
                    generatePgSA();              
                
//                this->lookupTable.generateFromSA(this, this->getElementsCount());
                buildReadsWithDuplicatesFilter();
            };

            SparseSuffixArray(PseudoGenome* pseudoGenome, std::istream& src)
            :   SparseSuffixArrayBase(determineElementsCount(pseudoGenome), getSizeInBytesWithGuard(pseudoGenome), pseudoGenome, pseudoGenome->getSymbolsPerElement()-1, sizeof(uint_skipped_element)),
                OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass, ThisDefaultSuffixArrayType>(
                            *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
                pseudoGenome(pseudoGenome),
                readsList(pseudoGenome->getReadsList()),
                sPacker(new SymbolsPackingFacility<uint_skipped_element>(pseudoGenome->getReadsSetProperties(), pseudoGenome->getSymbolsPerElement()-1)),
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
            
            ~SparseSuffixArray() {
                delete[] (this->suffixArray);
                delete(pseudoGenome);
                delete(sPacker);
            };
            
            inline OccurrencesIterator& getKmerOccurrencesIteratorImpl(const string& kmer) {
                this->findOccurrencesOf(kmer.data(), kmer.length());
                return *this;
            };
            
            char_pg kmerByPos[USHRT_MAX];
            
            inline OccurrencesIterator& getKmerOccurrencesIteratorImpl(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength) {
//                const char_pg* kmer = this->pseudoGenome->getSuffixPtrByPosition(originalIdx, pos);
                this->pseudoGenome->getKmerByPosition(originalIdx, pos, kmerLength, kmerByPos);
                this->findOccurrencesOf(kmerByPos, kmerLength);
                return *this;
            };
            
            inline void findKmerRangeImpl(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
                if (!lookupTable.findSARange(kmer, kmerLength, range))
                    kmerRangeBSearch(kmer, kmerLength, range);
                
                if ((checkLastElementFromShift > shift) && (range.end == this->getElementsCount()))
                    range.end--;   
            }
            
            const char* kmer;
            uint_read_len kmerLength;
            SARange<uint_pg_len> range;
            uchar shift;
            
            inline bool matchPrefix(const sa_pos_addr& saPosAddress) {
                for(int i = 0; i < shift; i++) 
                    if (kmer[i] != sPacker->reverseValue(getSkippedElementByAddress(saPosAddress), i + skippedSymbolsCount - shift))
                        return false;
                        
                return true;
            }
            
            inline void findOccurrencesOfImpl(const char* kmer, const uint_read_len& kmerLength) {
                this->kmer = kmer;
                this->kmerLength = kmerLength;
                this->shift = 0;
                determineCheckLastElementFromShift(this->kmerLength);
                this->findKmerRange(this->kmer, this->kmerLength, range);
                this->readsIterator.initIteration(kmerLength);
                if (range.start != range.end)
                    this->readsIterator.setIterationPosition(this->getPosition(range.start), 0);
            }
            
            inline void findOccurrences(uchar shift) { 
                this->shift = shift;
                const char* shiftedKmer = kmer + shift;
                uint_read_len shiftedLength = kmerLength - shift;
                this->findKmerRange(shiftedKmer, shiftedLength, range);

                sa_pos_addr posAddress = 0;
                while(range.start < range.end && !matchPrefix(posAddress = saPosIdx2Address(range.start)))
                    range.start++;
                if (range.start >= range.end) {
                    if (shift < skippedSymbolsCount)
                        findOccurrences(++shift);
                } else 
                    this->readsIterator.setIterationPosition(this->getPositionByAddress(posAddress), shift);
            };

            inline bool moveNextImpl() {
                if (this->readsIterator.moveNext())
                    return true;
                
                sa_pos_addr posAddress;
                while (true) {
                    while (++(range.start) < range.end) {
                        if(matchPrefix(posAddress = saPosIdx2Address(range.start))) {
                            this->readsIterator.setIterationPosition(this->getPositionByAddress(posAddress), shift);
                            if (this->readsIterator.moveNext())
                                return true;
                        }
                    }
                    if (shift == skippedSymbolsCount)
                        return false;
                    findOccurrences(++shift);
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
            
            string getTypeID() { return PGSATYPE_SPARSE; };

            inline const string getDescriptionImpl() { 
                string desc;
                desc = desc + "Sparse SA (interval = " + toString(this->getSkippedSymbolsCount()+1);
                uint_max size = sizeof(this) + this->getSizeInBytes() + 
                        sizeof(this->lookupTable) + this->lookupTable.getLookupTableLengthWithGuard() * sizeof(uint_pg_len) +
                        sizeof(this->pseudoGenome) + this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement() +
                        sizeof(this->readsList) + (uint_max) this->readsList->getReadsCount() * this->readsList->getListElementSize();
                desc = desc + ")\tTOTAL (Pg+RL+SA+LT) " + toMB(size, 2) + " MB\n" +
                        + "Pg: " + toMB(this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement(), 2) + " MB\t" 
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