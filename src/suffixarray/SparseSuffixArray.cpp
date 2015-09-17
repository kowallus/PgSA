#include "SparseSuffixArray.h"
#include "SparseSuffixArrayFactory.h"

#include "../sais/sais.h"

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::SparseSuffixArray(PseudoGenome* pseudoGenome, int fixed_min_k)
    : // up to 2 additional bytes to avoid overflowing during casting to uint_reads_cnt
    SparseSuffixArrayBase(determineElementsCount(pseudoGenome), getSuffixArraySizeInBytesWithGuard(pseudoGenome->getLength() / pseudoGenome->getSymbolsPerElement()), pseudoGenome, pseudoGenome->getSymbolsPerElement() - 1, sizeof (uint_skipped_element)),
    OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass, ThisDefaultSuffixArrayType>(
    *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
    pseudoGenome(pseudoGenome),
    readsList(pseudoGenome->getReadsList()),
    suffixArray(new uchar[this->getSizeInBytes()]),
    sPacker(new SymbolsPackingFacility<uint_skipped_element>(pseudoGenome->getReadsSetProperties(), pseudoGenome->getSymbolsPerElement() - 1)),
    lookupTable(lookupTableKeyPrefixLength, pseudoGenome->getReadsSetProperties()),
    fixed_min_k(fixed_min_k) {
        //TODO: check if setting 0 is necessary.... 
        suffixArray[this->getSizeInBytes() - 2] = 0;
        suffixArray[this->getSizeInBytes() - 1] = 0;

        pgStatic = this->pseudoGenome;
        maxReadLengthRaw = (this->pseudoGenome->maxReadLength() + this->pseudoGenome->getSymbolsPerElement() - 1) / this->pseudoGenome->getSymbolsPerElement();

        if ((sizeof (uint_pg_element) == sizeof (char)))
            generateSaisPgSA();
        else
            generatePgSA();

        this->lookupTable.generateFromPg(this->pseudoGenome, this->readsList, pseudoGenome->getSymbolsPerElement(), skippedSymbolsCount, fixed_min_k);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::prepareUnsortedSA() {

        const uchar* curSAPos = suffixArray;

        uint_reads_cnt readsListIndex = 0;

        for (uint_pg_len pgPos = 0; pgPos <= pseudoGenome->getLength(); pgPos += pseudoGenome->getSymbolsPerElement()) {
            uint_pg_len pgSuffixPos = pgPos + skippedSymbolsCount;
            while ((readsListIndex < readsList->getReadsCount()) &&
                    (readsList->getReadPosition(readsListIndex + 1) <= pgSuffixPos))
                readsListIndex++;
            *((uint_reads_cnt*) curSAPos) = readsListIndex;
            *((uint_read_len*) (curSAPos + POS_START_OFFSET)) = pgSuffixPos - readsList->getReadPosition(readsListIndex);
            *((uint_skipped_element*) (curSAPos + SKIPPED_OFFSET)) =
                    sPacker->packSymbols(pseudoGenome->getSuffix(pgPos, skippedSymbolsCount).data());

            curSAPos += SA_ELEMENT_SIZE;
        }

        elementsCountWithoutGuard = this->getPosition(this->elementsCount - 1).readListIndex == this->pseudoGenome->readsCount() ? this->elementsCount - 1 : this->elementsCount;

        if (curSAPos != suffixArray + this->elementsCount * (uint_max) SA_ELEMENT_SIZE)
            cout << "WARNING: SA generation failed: " << (int) (curSAPos - suffixArray) / SA_ELEMENT_SIZE << " elements instead of " << this->elementsCount << "\n";

    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::updateSAPositionQueue(int saGroup) {
        int relativePosition = 0;

        while (saPartSrc[saGroup]->read((char*) &relativePosition, sizeof (int))) {

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

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::generateSaisPgSA() {
        clock_checkpoint();

        int maxPartBruttoSize = INT_MAX / (sizeof (int));
        uint_pg_len pgElementsCountWithGuard = pseudoGenome->getElementsCountWithGuard();

        if ((uint_max) maxPartBruttoSize * 6 > pgElementsCountWithGuard * (uint_max) SA_ELEMENT_SIZE)
            maxPartBruttoSize = pgElementsCountWithGuard * (uint_max) SA_ELEMENT_SIZE / 6;

        maxPartSize = maxPartBruttoSize - maxReadLengthRaw;
        int noOfParts = 0;

        int* saisSA = (int*) malloc((size_t) (maxPartBruttoSize * sizeof (int)));

        while (noOfParts * (uint_pg_len) maxPartSize < pgElementsCountWithGuard) {

            uint_pg_len partSize = pgElementsCountWithGuard - noOfParts * (uint_pg_len) maxPartSize;
            if (partSize > (uint_pg_len) maxPartBruttoSize)
                partSize = maxPartBruttoSize;

            if (sais((const unsigned char*) this->pseudoGenome->getRawSuffix(1 + noOfParts * (uint_pg_len) maxPartSize), saisSA, (int) partSize) != 0) {
                fprintf(stderr, "Cannot allocate memory.\n");
                exit(EXIT_FAILURE);
            }

            std::ofstream dest("saisSA" + toString(noOfParts++) + ".tmp", std::ios::out | std::ios::binary);
            PgSAHelpers::writeArray(dest, saisSA, (size_t) (partSize) * sizeof (int));
            dest.close();
        }

        free((void*) saisSA);

        cout << "SAIS generation time " << clock_millis() << " msec!\n";

        clock_checkpoint();

        prepareUnsortedSA();

        readsList->buildLUT();

        currentSaPartPos.resize(noOfParts);
        for (int i = 0; i < noOfParts; i++) {
            saPartSrc.push_back(new std::ifstream("saisSA" + toString(i) + ".tmp", std::ifstream::binary));
            updateSAPositionQueue(i);
        }

        uint_read_len maxReadLength = this->pseudoGenome->maxReadLength();
        const uchar* curSAPos = suffixArray;
        uint_pg_len saPos;
        
        uint_pg_len i = 0;
       
        while (!saPartsOrder.empty()) {
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
                copyElementsByAddress((sa_pos_addr*) curSAPos, saPosIdx2Address(saPos));

            if (*(curSAPos + POS_START_OFFSET) <= maxReadLength - fixed_min_k + skippedSymbolsCount) {
                curSAPos += SA_ELEMENT_SIZE;
                i++;
            }

            saPartsOrder.pop_front();
            updateSAPositionQueue(j);
        }

        elementsCount = i;
        suffixArrayBytes = getSuffixArraySizeInBytesWithGuard(elementsCount);
        
        for (int i = 0; i < noOfParts; i++) {
            saPartSrc[i]->close();
            delete(saPartSrc[i]);
            remove(("saisSA" + toString(noOfParts) + ".tmp").c_str());
        }

        saPartSrc.clear();

        cout << "SA generation time " << clock_millis() << " msec!\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::swapElementsByAddress(const sa_pos_addr saPosAddressFst, const sa_pos_addr saPosAddressSnd) {
        uchar tmp;
        for (int i = 0; i < SA_ELEMENT_SIZE; i++) {
            tmp = *(((uchar*) saPosAddressFst) + i);
            *(((uchar*) saPosAddressFst) + i) = *(((uchar*) saPosAddressSnd) + i);
            *(((uchar*) saPosAddressSnd) + i) = tmp;
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::copyElementsByAddress(const sa_pos_addr saPosAddressDest, const sa_pos_addr saPosAddressSrc) {
        std::copy((uchar*) saPosAddressSrc, (uchar*) saPosAddressSrc + SA_ELEMENT_SIZE, (uchar*) saPosAddressDest);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::generatePgSA() {
        clock_checkpoint();

        prepareUnsortedSA();

        qsort(suffixArray, elementsCountWithoutGuard, sizeof (uchar) * SA_ELEMENT_SIZE, this->pgSuffixesCompare);

        if (fixed_min_k > 1) {
            const uchar* curSAPos = suffixArray;
            const uchar* oldSAPos = suffixArray;
            uint_read_len maxReadLength = this->pseudoGenome->maxReadLength();
            
            for(uint_pg_len i = 0; i < elementsCountWithoutGuard; i++) {
                if (*(oldSAPos + POS_START_OFFSET) <= maxReadLength - fixed_min_k + skippedSymbolsCount) {
                    copyElementsByAddress((sa_pos_addr*) curSAPos, (sa_pos_addr*) oldSAPos);
                    curSAPos += SA_ELEMENT_SIZE;
                }    
                oldSAPos += SA_ELEMENT_SIZE;
            }
            
            elementsCount = (curSAPos - suffixArray) / SA_ELEMENT_SIZE;
            suffixArrayBytes = getSuffixArraySizeInBytesWithGuard(elementsCount);
        }
        
        cout << "SA generation time " << clock_millis() << " msec!\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline int SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::pgSuffixesCompare(const void* a, const void* b) {
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
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::SparseSuffixArray(PseudoGenome* pseudoGenome, std::istream& src)
    : SparseSuffixArrayBase(determineElementsCount(pseudoGenome), getSuffixArraySizeInBytesWithGuard(pseudoGenome->getLength() / pseudoGenome->getSymbolsPerElement()), pseudoGenome, pseudoGenome->getSymbolsPerElement() - 1, sizeof (uint_skipped_element)),
    OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename ReadsListIteratorFactoryTemplate<ReadsListClass>::ReadsListIteratorClass, ThisDefaultSuffixArrayType>(
    *ReadsListIteratorFactoryTemplate<ReadsListClass>::getReadsListIterator(* ((ReadsListClass*) pseudoGenome->getReadsList()))),
    pseudoGenome(pseudoGenome),
    readsList(pseudoGenome->getReadsList()),
    sPacker(new SymbolsPackingFacility<uint_skipped_element>(pseudoGenome->getReadsSetProperties(), pseudoGenome->getSymbolsPerElement() - 1)),
    lookupTable(lookupTableKeyPrefixLength, pseudoGenome->getReadsSetProperties()) {
        uint_max arraySize;
        src >> arraySize;
        src.get(); // '/n'

        if (arraySize > this->getSizeInBytes())
            cout << "WARNING: wrong size of suffixarray.";

        suffixArrayBytes = arraySize;
        elementsCount = getSuffixArrayElementsCount(suffixArrayBytes);
        
        suffixArray = (uchar*) PgSAHelpers::readArray(src, this->getSizeInBytes() * sizeof (uchar));
        src.get(); // '/n'

        this->lookupTable.read(src);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::~SparseSuffixArray() {
        delete[] (this->suffixArray);
        delete(pseudoGenome);
        delete(sPacker);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::write(ostream& dest) {
        dest << this->getSizeInBytes() << "\n";
        PgSAHelpers::writeArray(dest, suffixArray, this->getSizeInBytes());
        dest << "\n";
        lookupTable.write(dest);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::findKmerRangeImpl(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) {
        if (!lookupTable.findSARange(kmer, kmerLength, range))
            kmerRangeBSearch(kmer, kmerLength, range);

        if ((checkLastElementFromShift > shift) && (range.end == this->getElementsCount()))
            range.end--;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline typename SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>
    ::int_kmers_comp SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::kmerSAComp(const uint_pg_element* kmerPtr, const uint_read_len& directCompareElementsCount, const uint_read_len& omitElementsCount, const uchar& symbolsLeftToCompare, const uint_pg_len& suffixIdx) {
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

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::kmerRangeBSearch(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) {
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
                if (range.start < mIdx)
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
                if (range.start < lIdx)
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

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::determineCheckLastElementFromShift(uint_read_len& kmerLength) {
        RPGOffset<uint_read_len, uint_reads_cnt> lastElement = this->getPosition(this->getElementsCount() - 1);
        if (lastElement.readListIndex == pseudoGenome->readsCount())
            checkLastElementFromShift = lastElement.offset + kmerLength;
        else
            checkLastElementFromShift = 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::findOccurrences(uchar shift) {
        this->shift = shift;
        const char* shiftedKmer = kmer + shift;
        uint_read_len shiftedLength = kmerLength - shift;
        this->findKmerRange(shiftedKmer, shiftedLength, range);

        sa_pos_addr posAddress = 0;
        while (range.start < range.end && !matchPrefix(posAddress = saPosIdx2Address(range.start)))
            range.start++;
        if (range.start >= range.end) {
            if (shift < skippedSymbolsCount)
                findOccurrences(++shift);
        } else
            this->readsIterator.setIterationPosition(this->getPositionByAddress(posAddress), shift);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline void SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::findOccurrencesOfImpl(const char* kmer, const uint_read_len& kmerLength) {
        this->kmer = kmer;
        this->kmerLength = kmerLength;
        this->shift = 0;
        determineCheckLastElementFromShift(this->kmerLength);
        this->findKmerRange(this->kmer, this->kmerLength, range);
        this->readsIterator.initIteration(kmerLength);
        if (range.start != range.end)
            this->readsIterator.setIterationPosition(this->getPosition(range.start), 0);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline typename SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>
    ::OccurrencesIterator& SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getKmerOccurrencesIteratorImpl(const string& kmer) {
        this->findOccurrencesOf(kmer.data(), kmer.length());
        return *this;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline string SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getKmerImpl(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength) {
        this->pseudoGenome->getKmerByPosition(originalIdx, pos, kmerLength, kmerByPos);
        return kmerByPos;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline bool SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::moveNextImpl() {
        if (this->readsIterator.moveNext())
            return true;

        sa_pos_addr posAddress;
        while (true) {
            while (++(range.start) < range.end) {
                if (matchPrefix(posAddress = saPosIdx2Address(range.start))) {
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
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline bool SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::matchPrefix(const sa_pos_addr& saPosAddress) {
        for (int i = 0; i < shift; i++)
            if (kmer[i] != sPacker->reverseValue(getSkippedElementByAddress(saPosAddress), i + skippedSymbolsCount - shift))
                return false;

        return true;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    const uint_max SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::determineElementsCount(PseudoGenome* pseudoGenome) {
        return pseudoGenome->getLength() / pseudoGenome->getSymbolsPerElement() + 1;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    uint_max SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getSuffixArraySizeInBytesWithGuard(uint_pg_len elementsCount) {
        return sizeof (uchar) * (elementsCount + 2) * SA_ELEMENT_SIZE;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    uint_pg_len SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getSuffixArrayElementsCount(uint_max sizeInBytesWithGuard) {
        return sizeInBytesWithGuard / sizeof(uchar) / SA_ELEMENT_SIZE - 2;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_read_len SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getPosStartOffsetByAddress(const sa_pos_addr saPosAddress) {
        return *((uint_read_len*) ((uchar*) saPosAddress + POS_START_OFFSET));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const RPGOffset<uint_read_len, uint_reads_cnt> SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getPositionByAddress(const sa_pos_addr& saPosAddress) {
        return
        {
            getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress)
        };
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_pg_element* SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getRawSuffixByAddress(const sa_pos_addr& saPosAddress) {
        return this->pseudoGenome->getRawSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_reads_cnt SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getReadsListIndexByAddress(const sa_pos_addr saPosAddress) {
        return (*((uint_reads_cnt*) (saPosAddress))) & READSLIST_INDEX_MASK;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_skipped_element SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getSkippedElementByAddress(const sa_pos_addr saPosAddress) {
        return *((uint_skipped_element*) ((uchar*) saPosAddress + SKIPPED_OFFSET));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const string SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getSuffixByAddress(const sa_pos_addr& saPosAddress, const uint_pg_len length) {
        return this->pseudoGenome->getSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress), length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const sa_pos_addr SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::saPosIdx2Address(const uint_pg_len posIdx) {
        return (sa_pos_addr) (suffixArray + ((uint_max) posIdx) * SA_ELEMENT_SIZE);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const RPGOffset<uint_read_len, uint_reads_cnt> SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getPositionImpl(const uint_pg_len posIdx) {
        return getPositionByAddress(saPosIdx2Address(posIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const string SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getSuffixImpl(const uint_pg_len posIdx, const uint_pg_len length) {
        return getSuffixByAddress(saPosIdx2Address(posIdx), length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline ReadsSetInterface<uint_read_len, uint_reads_cnt>* SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getReadsSetImpl() {
        return this->pseudoGenome;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    const string SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getDescriptionImpl() {
        string desc;
        desc = desc + "Sparse SA (interval = " + toString(this->getSkippedSymbolsCount() + 1);
        uint_max size = sizeof (this) + this->getSizeInBytes() +
                sizeof (this->lookupTable) + this->lookupTable.getLookupTableLengthWithGuard() * sizeof (uint_pg_len) +
                sizeof (this->pseudoGenome) + this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement() +
                sizeof (this->readsList) + (uint_max) this->readsList->getReadsCount() * (this->readsList->getListElementSize() + sizeof(uint_reads_cnt));
        desc = desc + ")\tTOTAL (Pg+RL+SA+LT) " + toMB(size, 2) + " MB\n" +
                +"Pg: " + toMB(this->pseudoGenome->getElementsCountWithGuard() * this->pseudoGenome->getBytesPerElement(), 2) + " MB\t"
                + "RL: " + toMB((uint_max) this->readsList->getReadsCount() * (this->readsList->getListElementSize() + sizeof(uint_reads_cnt)), 2) + " MB\t"
                + "SA: " + toMB(this->getSizeInBytes(), 2) + " MB\t"
                + "LT: " + toMB(this->lookupTable.getLookupTableLengthWithGuard() * sizeof (uint_pg_len), 2) + " MB\n";
        desc = desc + toString(this->readsList->getReadsCount()) + " ";
        if (this->pseudoGenome->isReadLengthConstant())
            desc = desc + "constant";
        else
            desc = desc + "variable";
        desc = desc + " length reads\t max length: " + toString(this->pseudoGenome->maxReadLength())
                + "\t Pg length: " + toString(this->pseudoGenome->getLength()) + "\n";
        return desc;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_pg_element* SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getRawSuffix(const uint_pg_len posIdx) {
        return getRawSuffixByAddress(saPosIdx2Address(posIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    inline const uint_pg_element* SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getRawSuffixStatic(const sa_pos_addr& saPosAddress) {
        return pgStatic->getRawSuffix(getReadsListIndexByAddress(saPosAddress), getPosStartOffsetByAddress(saPosAddress));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, uchar SA_ELEMENT_SIZE, uchar POS_START_OFFSET, uchar SKIPPED_OFFSET, uint_reads_cnt READSLIST_INDEX_MASK, class ReadsListClass>
    string SparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, SA_ELEMENT_SIZE, POS_START_OFFSET, SKIPPED_OFFSET, READSLIST_INDEX_MASK, ReadsListClass>::getTypeID() {
        return PGSATYPE_SPARSE;
    }
    
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    
    // unused variants
    // TODO: consider resigning from uint_ps_element_std for skipped elements or even maybe for Pg...

    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class SparseSuffixArray<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;

    
}
