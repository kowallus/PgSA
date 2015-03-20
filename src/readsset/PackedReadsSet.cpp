#include "PackedReadsSet.h"

namespace PgSAReadsSet {
        
    template<class ReadsSourceIterator>
    PackedReadsSet::PackedReadsSet(ReadsSourceIterator* readsIterator) {

        bool symbolOccured[UCHAR_MAX] = {0};

        while (readsIterator->moveNext()) {

            properties->readsCount++;

            // analize read length
            uint_read_len_max length = readsIterator->getReadLength();
            if (properties->maxReadLength == 0)
                properties->maxReadLength = length;
            else if (properties->maxReadLength != length) {
                properties->constantReadLength = false;
                if (properties->maxReadLength < length)
                    properties->maxReadLength = length;
            }

            properties->allReadsLength += length;

            //analize symbols
            string read(readsIterator->getRead());
            
            for (uint_read_len_max i = 0; i < length; i++) {
                read[i] = toupper(read[i]);
                if (!symbolOccured[(unsigned char) read[i]]) {
                    symbolOccured[(unsigned char) read[i]] = true;
                    properties->symbolsCount++;
                }
            }
            
            lengths.push_back((uint_read_len_max) read.length());
        }

        // order symbols

        for (uint_symbols_cnt i = 0, j = 0; i < properties->symbolsCount; j++)
            if (symbolOccured[(unsigned char) j])
                properties->symbolsList[(unsigned char) (i++)] = j;

        properties->generateSymbolOrder();

        sPacker = new SymbolsPackingFacility<uint_ps_element_min>(properties, SymbolsPackingFacility<uint_ps_element_min>::maxSymbolsPerElement(properties->symbolsCount));
        unsigned char buffer[properties->maxReadLength + 1] = {0}; 

        readsIterator->rewind();
        while (readsIterator->moveNext()) {
            int packedLength = sPacker->packSequence(readsIterator->getRead().c_str(), readsIterator->getReadLength(), buffer);
            packedReads.push_back(vector<uint_ps_element_min>(buffer, buffer + packedLength));
        }
               
    };

    int PackedReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx){
        return sPacker->compareSequences(packedReads[lIdx].data(), packedReads[rIdx].data(), properties->maxReadLength);
    }

    int PackedReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx, uint_read_len_max offset) {
        return sPacker->compareSequences(packedReads[lIdx].data(), packedReads[rIdx].data(), offset, properties->maxReadLength - offset);
    }

    int PackedReadsSet::compareSuffixWithPrefix(uint_reads_cnt_max sufIdx, uint_reads_cnt_max preIdx, uint_read_len_max sufOffset) {
        return sPacker->compareSuffixWithPrefix(packedReads[sufIdx].data(), packedReads[preIdx].data(), sufOffset, properties->maxReadLength - sufOffset);
    }

    template PackedReadsSet::PackedReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_min>>(ConcatenatedReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_std>>(ConcatenatedReadsSourceIterator<uint_read_len_std>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTAReadsSourceIterator<uint_read_len_min>>(FASTAReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTAReadsSourceIterator<uint_read_len_std>>(FASTAReadsSourceIterator<uint_read_len_std>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTQReadsSourceIterator<uint_read_len_min>>(FASTQReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTQReadsSourceIterator<uint_read_len_std>>(FASTQReadsSourceIterator<uint_read_len_std>* readsIterator);
    
}
