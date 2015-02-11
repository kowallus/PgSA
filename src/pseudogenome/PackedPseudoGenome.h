#ifndef PACKEDPSEUDOGENOME_H_INCLUDED
#define PACKEDPSEUDOGENOME_H_INCLUDED

#include "DefaultPseudoGenome.h"
#include "PackedPseudoGenomeBase.h"
#include "readslist/ReadsListTypes.h"
#include "packing/SymbolsPackingFacility.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    const string PGTYPE_PACKED = "PACKED_PGEN";

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass >
    class PackedPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len>,
            public PackedPseudoGenomeBase
    {
        private:
            uint_pg_element* sequence;
            const char_pg* orgPg;
            
            SymbolsPackingFacility<uint_pg_element>* sPacker;
            
            PseudoGenomeBase* srcPGB = 0; //if set manages readsList;
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

            CountQueriesCacheBase* countQueriesCache = 0;

            inline const uint_pg_element* getRawSuffix(const uint_pg_len rawPos) { return sequence + rawPos; };
            
        public:

            PackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement)
            : PackedPseudoGenomeBase(srcPseudoGenome->getLength(), srcPseudoGenome->getReadsSetProperties(), symbolsPerElement, sizeof(uint_pg_element))
            {
                this->sequence = new uint_pg_element[this->getElementsCountWithGuard()];
                sPacker = new SymbolsPackingFacility<uint_pg_element>(this->getReadsSetProperties(), symbolsPerElement);
                
                // first element is shifted one symbol...
                sequence[0] = sPacker->packPrefixSymbols(srcPseudoGenome->getSuffix(0), symbolsPerElement - 1);
                
                uint_max i = sPacker->packSequence(srcPseudoGenome->getSuffix(symbolsPerElement - 1), srcPseudoGenome->getLengthWithGuard(), sequence + 1);
                while (++i < getElementsCountWithGuard())
                    sequence[i] = sPacker->getMaxValue();
                
                this->readsList = srcPseudoGenome->getReadsList();
                
                countQueriesCache = srcPseudoGenome->getCountQueriesCacheBase();
                srcPGB = srcPseudoGenome;
                orgPg = srcPseudoGenome->getSuffix(0);
            }

            PackedPseudoGenome(uint_pg_len pgLength, std::istream& src)
            : PackedPseudoGenomeBase(pgLength, src, sizeof(uint_pg_element)){
                
                uint_pg_len arraySize;
                src >> arraySize;
                src.get(); // '/n'

                if (arraySize != getElementsCountWithGuard())
                    cout << "WARNING: wrong size of pseudogenome.";
                
                sequence = (uint_pg_element*) PgSAHelpers::readArray(src, getElementsCountWithGuard() * sizeof(uint_pg_element));
                src.get(); // '/n'
                this->readsList = new ReadsListClass(maxReadLength(), src);
                
                sPacker = new SymbolsPackingFacility<uint_pg_element>(this->getReadsSetProperties(), symbolsPerElement);
                
                string str = getSuffix(0, this->getLength());
                char* seq = new char_pg[str.size() + 1];
                std::copy(str.begin(), str.end(), seq);
                seq[str.size()] = '\0';
                orgPg = seq;
            }

            ~PackedPseudoGenome() {
                delete[]sequence;
                if (srcPGB == 0)
                    delete(readsList);
                else
                    delete[]orgPg;
                delete(sPacker);
            }

            static PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* castBase(PseudoGenomeBase* base) {
                // TODO: validate
                return static_cast<PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>*>(base);
            }

            string getTypeID() { return PGTYPE_PACKED; };

            void write(std::ostream& dest) {
                //FIXME: delegate to base type
                this->getReadsSetProperties()->write(dest);
                dest << (int) this->symbolsPerElement << "\n";
                
                dest << getElementsCountWithGuard() << "\n";
                PgSAHelpers::writeArray(dest, sequence, getElementsCountWithGuard() * sizeof(uint_pg_element));
                dest << "\n";
                this->readsList->write(dest);
            }
            
            CountQueriesCacheBase* getCountQueriesCacheBase() {
                return countQueriesCache;
            }
            
            uint_pg_len getElementsCountWithGuard() { return 2 + (this->length + this->properties->maxReadLength - 1) / symbolsPerElement; };
            
            inline ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList() { return readsList; }

//            inline uint_pg_len getLength() { return this->length; };

            inline char getSymbol(uint_pg_len pos) { 
                uint_max division = divideBySmallInteger(pos + 1, symbolsPerElement);
                return sPacker->reverseValue(sequence[division], moduloBySmallInteger(pos + 1, symbolsPerElement, division)); 
            };

            inline const uint_pg_element* getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos) { return getRawSuffix((this->readsList->getReadPosition(readsListIdx) + pos + 1) / symbolsPerElement); };
            
            inline uint_read_len maxReadLength() { return this->properties->maxReadLength; };

            inline uint_reads_cnt readsCount() { return this->properties->readsCount; };

            inline bool isReadLengthConstant() { return this->properties->constantReadLength; };

            inline const string getSuffix(const uint_pg_len pos, const uint_pg_len length) { 
                return sPacker->reverseSequence(sequence, pos + 1, length);
            };

            inline const string getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length) { 
                return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset, length); 
            };
            
            inline const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
                return orgPg + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
            }
            
            inline void getKmerByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kValue, char_pg* kmerPtr) {
                sPacker->reverseSequence(sequence, this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos + 1, kValue, kmerPtr);
            }
            
            const string getRead(uint_reads_cnt originalIdx) {               
                return getSuffix(this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)), readLength(originalIdx));
            };
            
            uint_read_len readLength(uint_reads_cnt originalIdx) { 
                return this->readsList->getReadLength(
                        this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)); 
            };

            SymbolsPackingFacility<uint_pg_element>* getSymbolsPacker() {
                return sPacker;
            }
            
            uint_read_len maxReadLengthVirtual() { return maxReadLength(); };
            uint_reads_cnt readsCountVirtual() { return readsCount(); };
            bool isReadLengthConstantVirtual() { return isReadLengthConstant(); };
            const string getReadVirtual(uint_reads_cnt i) { return getRead(i); };
            uint_read_len readLengthVirtual(uint_reads_cnt i) { return readLength(i); };
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    using PackedPseudoGenomeOfConstantLengthReadsType = PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
};

#endif // PACKEDPSEUDOGENOME_H_INCLUDED
