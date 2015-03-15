#include "MultiPackedPseudoGenome.h"

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::MultiPackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement)
    : PackedPseudoGenomeBase(srcPseudoGenome->getLength(), srcPseudoGenome->getReadsSetProperties(), symbolsPerElement, sizeof (uint_pg_element)) {
        this->sequences = new uint_pg_element*[symbolsPerElement];
        sPacker = new SymbolsPackingFacility<uint_pg_element>(this->getReadsSetProperties(), symbolsPerElement);

        for (int j = 0; j < symbolsPerElement; j++) {
            this->sequences[j] = new uint_pg_element[this->getElementsCountWithGuard()];

            uint_max i = sPacker->packSequence(srcPseudoGenome->getSuffix(j), srcPseudoGenome->getLengthWithGuard(), sequences[j]);
            while (++i < getElementsCountWithGuard())
                sequences[j][i] = sPacker->getMaxValue();
        }

        this->readsList = srcPseudoGenome->getReadsList();
        countQueriesCache = srcPseudoGenome->getCountQueriesCacheBase();
        srcPGB = srcPseudoGenome;
        orgPg = srcPseudoGenome->getSuffix(0);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::~MultiPackedPseudoGenome() {
        for (int i = 0; i < symbolsPerElement; i++)
            delete[]sequences[i];
        delete[]sequences;
        delete(sPacker);
        //delete(readsList);
        delete(srcPGB);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    void MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::write(std::ostream& dest) {
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const uint_pg_element* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRawSuffix(const uint_pg_len pos) {
        uint_pg_len division = divideBySmallInteger(pos, symbolsPerElement);
        return sequences[moduloBySmallInteger(pos, symbolsPerElement, division)] + division;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const uint_pg_element* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos) {
        return getRawSuffix(this->readsList->getReadPosition(readsListIdx) + pos);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSuffix(const uint_pg_len pos, const uint_pg_len length) {
        uint_pg_len division = divideBySmallInteger(pos, symbolsPerElement);
        uint_pg_len reminder = moduloBySmallInteger(pos, symbolsPerElement, division);
        return sPacker->reverseSequence(sequences[reminder], division * symbolsPerElement, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length) {
        return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const char_pg* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
        return orgPg + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    char MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSymbol(uint_pg_len pos) {
        uint_max division = divideBySmallInteger(pos, symbolsPerElement);
        return sPacker->reverseValue(sequences[0][division], moduloBySmallInteger((uint_max) pos, symbolsPerElement, division));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getReadsList() {
        return readsList;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRead(uint_reads_cnt originalIdx) {
        return getSuffix(this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)), readLength(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_pg_len MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getElementsCountWithGuard() {
        return 2 + (this->length + this->properties->maxReadLength - 1) / symbolsPerElement;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    string MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getTypeID() {
        return PGTYPE_MULTIPACKED;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    CountQueriesCacheBase* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getCountQueriesCacheBase() {
        return countQueriesCache;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    SymbolsPackingFacility<uint_pg_element>* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSymbolsPacker() {
        return sPacker;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::castBase(PseudoGenomeBase* base) {
        // TODO: validate
        return static_cast<MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>*> (base);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    bool MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::isReadLengthConstant() {
        return this->properties->constantReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::maxReadLength() {
        return this->properties->maxReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readLength(uint_reads_cnt originalIdx) {
        return this->readsList->getReadLength(
                this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_reads_cnt MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readsCount() {
        return this->properties->readsCount;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getReadVirtual(uint_reads_cnt i) {
        return getRead(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    bool MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::isReadLengthConstantVirtual() {
        return isReadLengthConstant();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::maxReadLengthVirtual() {
        return maxReadLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readLengthVirtual(uint_reads_cnt i) {
        return readLength(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_reads_cnt MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readsCountVirtual() {
        return readsCount();
    }
    
    template class MultiPackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;

    template class MultiPackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class MultiPackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    
}
