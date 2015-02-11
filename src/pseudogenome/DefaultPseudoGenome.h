#ifndef DEFAULTPSEUDOGENOME_H_INCLUDED
#define DEFAULTPSEUDOGENOME_H_INCLUDED

#include "PseudoGenomeInterface.h"
#include "readslist/ReadsListTypes.h"
#include "../index/cache/CountQueriesCacheBase.h"
#include "PseudoGenomeBase.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

namespace PgSAIndex {

    const string PGTYPE_DEFAULT = "DEFAULT_PGEN";

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    class DefaultPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len>,
                               public PseudoGenomeBase
    {
        protected:
            char_pg* sequence;

            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

            CountQueriesCacheBase* countQueriesCache = 0;
            
        public:

            DefaultPseudoGenome(uint_pg_len pgLength, ReadsSetProperties* properties)
            : PseudoGenomeBase(pgLength, properties) {
                this->sequence = new char_pg[(this->getLengthWithGuard()) * sizeof(char_pg)];
            }

            DefaultPseudoGenome(uint_pg_len pgLength, std::istream& src)
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

            ~DefaultPseudoGenome() {
                delete[]sequence;
                delete(readsList);
            }

            static DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* castBase(PseudoGenomeBase* base) {
                // TODO: validate
                return static_cast<DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>*>(base);
            }

            string getTypeID() { return PGTYPE_DEFAULT; };

            void write(std::ostream& dest) {
                this->getReadsSetProperties()->write(dest);
                dest << getLengthWithGuard() << "\n";
                PgSAHelpers::writeArray(dest, sequence, getLengthWithGuard() * sizeof(char_pg));
                dest << "\n";
                this->readsList->write(dest);
            }
            
            CountQueriesCacheBase* getCountQueriesCacheBase() {
                return countQueriesCache;
            }
            
            inline ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList() { return readsList; }

//            inline uint_pg_len getLength() { return this->length; };

            uint_pg_len getLengthWithGuard() { return this->length + this->properties->maxReadLength; };
            
            inline char getSymbol(uint_pg_len pos) { return sequence[pos]; };

            inline const char_pg* getSuffix(const uint_pg_len pos) { return sequence + pos; };
            
            inline const char_pg* getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset) { return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset); };

            inline const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
                return sequence + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
            }
            
            inline uint_read_len maxReadLength() { return this->properties->maxReadLength; };

            inline uint_reads_cnt readsCount() { return this->properties->readsCount; };

            inline bool isReadLengthConstant() { return this->properties->constantReadLength; };

            const string getRead(const uint_reads_cnt originalIdx) {
                uint_pg_len pos = this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
                return string(sequence + pos, readLength(originalIdx));
            };
            
            uint_read_len readLength(const uint_reads_cnt originalIdx) { return this->readsList->getReadLength(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)); };

            uint_read_len maxReadLengthVirtual() { return maxReadLength(); };
            uint_reads_cnt readsCountVirtual() { return readsCount(); };
            bool isReadLengthConstantVirtual() { return isReadLengthConstant(); };
            const string getReadVirtual(uint_reads_cnt i) { return getRead(i); };
            uint_read_len readLengthVirtual(uint_reads_cnt i) { return readLength(i); };

    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using DefaultPseudoGenomeOfConstantLengthReadsType = DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass >
    class GeneratedPseudoGenome: public DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass >
    {
        private:
            uint_pg_len pos = 0;
            GeneratedReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>* genReadsList = 0;

        public:

            GeneratedPseudoGenome(uint_pg_len sequenceLength, ReadsSetProperties* properties)
            : DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass >(sequenceLength, properties)
            {
                this->pos = 0;

                this->genReadsList = new GeneratedReadsListClass (this->properties->maxReadLength, this->properties->readsCount, this->length);

                this->readsList = genReadsList;
            }

            ~GeneratedPseudoGenome() { }

            void append(const string& read, uint_read_len length, uint_read_len overlap, uint_reads_cnt orgIdx) {

                genReadsList->add(pos, length, orgIdx);

                uint_read_len len = length - overlap;
                if (len > 0) {
                    strncpy(this->sequence + pos, read.data(), len);
                    pos += len;
                }
            }

            void validate() {

                if (pos != this->length)
                    cout << "WARNING: Generated " << pos << " pseudogenome instead of " << this->length << "\n";

                genReadsList->validate();

                // adding guard
                for(uint_pg_len i = pos; i < this->getLengthWithGuard(); i++)
                    this->sequence[i] = this->properties->symbolsList[this->properties->symbolsCount - 1];
            }
            
            void setCountQueriesCache(CountQueriesCacheBase* cqcb) {
                this->countQueriesCache = cqcb;
            }
            
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using GeneratedPseudoGenomeOfConstantLengthReadsType = GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

};

#endif // DEFAULTPSEUDOGENOME_H_INCLUDED
