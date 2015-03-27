#include "DefaultPseudoGenome.h"

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultPseudoGenome(uint_pg_len pgLength, ReadsSetProperties* properties)
    : PseudoGenomeBase(pgLength, properties) {
        this->sequence = new char_pg[this->getLengthWithGuard()];
    }
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultPseudoGenome(uint_pg_len pgLength, std::istream& src)
    : PseudoGenomeBase(pgLength, src) {

        uint_pg_len arraySize;
        src >> arraySize;
        src.get(); // '/n'

        if (arraySize != getLengthWithGuard())
            cout << "WARNING: wrong size of pseudogenome.";

        sequence = (char_pg*) PgSAHelpers::readArray(src, getLengthWithGuard() * sizeof(char_pg));
        src.get(); // '/n'
        this->readsList = new ReadsListClass(maxReadLength(), src);
    }

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::~DefaultPseudoGenome() {
        if (readsList)
            delete(readsList);
        delete[]sequence;      
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::castBase(PseudoGenomeBase* base) {
        // TODO: validate
        return static_cast<DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>*> (base);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getTypeID() {
        return PGTYPE_DEFAULT;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::write(std::ostream& dest) {
        this->getReadsSetProperties()->write(dest);
        dest << getLengthWithGuard() << "\n";
        PgSAHelpers::writeArray(dest, sequence, getLengthWithGuard() * sizeof(char_pg));
        dest << "\n";
        this->readsList->write(dest);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getReadsList() 
    {
        return readsList; 
    }
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_pg_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getLengthWithGuard() {
        return this->length + this->properties->maxReadLength;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    char DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSymbol(uint_pg_len pos) {
        return sequence[pos];
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffix(const uint_pg_len pos) {
        return sequence + pos;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset) {
        return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
        return sequence + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::maxReadLength() {
        return this->properties->maxReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_reads_cnt DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readsCount() {
        return this->properties->readsCount;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isReadLengthConstant() {
        return this->properties->constantReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getRead(const uint_reads_cnt originalIdx) {
        uint_pg_len pos = this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
        return string(sequence + pos, readLength(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readLength(const uint_reads_cnt originalIdx) {
        return this->readsList->getReadLength(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const char_pg DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSymbolImpl(const uint_pg_len posIdx) {
        return *(sequence + posIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const uint_pg_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getLengthImpl() {
        return this->getPseudoGenomeLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getReadVirtual(uint_reads_cnt i) {
        return getRead(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readLengthVirtual(uint_reads_cnt i) {
        return readLength(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isReadLengthConstantVirtual() {
        return isReadLengthConstant();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_reads_cnt DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readsCountVirtual() {
        return readsCount();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::maxReadLengthVirtual() {
        return maxReadLength();
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::GeneratedPseudoGenome(uint_pg_len sequenceLength, ReadsSetProperties* properties)
    : DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass >(sequenceLength, properties) {
        this->pos = 0;

        this->genReadsList = new GeneratedReadsListClass(this->properties->maxReadLength, this->properties->readsCount, this->length);

        this->readsList = genReadsList;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::~GeneratedPseudoGenome() {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::append(const string& read, uint_read_len length, uint_read_len overlap, uint_reads_cnt orgIdx) {

        genReadsList->add(pos, length, orgIdx);

        uint_read_len len = length - overlap;
        if (len > 0) {
            strncpy(this->sequence + pos, read.data(), len);
            pos += len;
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::validate() {

        if (pos != this->length)
            cout << "WARNING: Generated " << pos << " pseudogenome instead of " << this->length << "\n";

        genReadsList->validate();

        // adding guard
        for (uint_pg_len i = pos; i < this->getLengthWithGuard(); i++)
            this->sequence[i] = this->properties->symbolsList[this->properties->symbolsCount - 1];
    }
    
    template class DefaultPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class DefaultPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    
    template class GeneratedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}
