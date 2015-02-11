#ifndef MULTIPACKEDPSEUDOGENOME_H_INCLUDED
#define MULTIPACKEDPSEUDOGENOME_H_INCLUDED

#include "DefaultPseudoGenome.h"
#include "PackedPseudoGenomeBase.h"
#include "readslist/ReadsListTypes.h"
#include "packing/SymbolsPackingFacility.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    const string PGTYPE_MULTIPACKED = "MULTIPACKED_PGEN";
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass >
    class MultiPackedPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len>,
            public PackedPseudoGenomeBase
    {
        private:
            uint_pg_element** sequences;
            const char_pg* orgPg;
            
            SymbolsPackingFacility<uint_pg_element>* sPacker;
            
            PseudoGenomeBase* srcPGB; //manages readsList;
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

            CountQueriesCacheBase* countQueriesCache = 0;
            
            inline const uint_pg_element* getRawSuffix(const uint_pg_len pos) {
                uint_pg_len division = divideBySmallInteger(pos, symbolsPerElement);
                return sequences[moduloBySmallInteger(pos, symbolsPerElement, division)] + division; 
            };
            
        public:

            MultiPackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement)
            : PackedPseudoGenomeBase(srcPseudoGenome->getLength(), srcPseudoGenome->getReadsSetProperties(), symbolsPerElement, sizeof(uint_pg_element))
            {
                this->sequences = new uint_pg_element*[symbolsPerElement];
                sPacker = new SymbolsPackingFacility<uint_pg_element>(this->getReadsSetProperties(), symbolsPerElement);
                
                for(int j = 0; j < symbolsPerElement; j++) {
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

            ~MultiPackedPseudoGenome() {
                for(int i = 0; i < symbolsPerElement; i++)
                    delete[]sequences[i];
                delete[]sequences;
                delete(sPacker);
                //delete(readsList);
                delete(srcPGB);
            }

            static MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* castBase(PseudoGenomeBase* base) {
                // TODO: validate
                return static_cast<MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>*>(base);
            }

            string getTypeID() { return PGTYPE_MULTIPACKED; };

            void write(std::ostream& dest) { }
            
            CountQueriesCacheBase* getCountQueriesCacheBase() {
                return countQueriesCache;
            }
            
            uint_pg_len getElementsCountWithGuard() { return 2 + (this->length + this->properties->maxReadLength - 1) / symbolsPerElement; };
            
            inline ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList() { return readsList; }

//            inline uint_pg_len getLength() { return this->length; };

            inline char getSymbol(uint_pg_len pos) { 
                uint_max division = divideBySmallInteger(pos, symbolsPerElement);
                return sPacker->reverseValue(sequences[0][division], moduloBySmallInteger(pos, symbolsPerElement, division)); 
            };

            // calculated position should be divisible by symbolsPerElement
            inline const uint_pg_element* getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos) { return getRawSuffix(this->readsList->getReadPosition(readsListIdx) + pos); };
            
            inline uint_read_len maxReadLength() { return this->properties->maxReadLength; };

            inline uint_reads_cnt readsCount() { return this->properties->readsCount; };

            inline bool isReadLengthConstant() { return this->properties->constantReadLength; };

            inline const string getSuffix(const uint_pg_len pos, const uint_pg_len length) {
                uint_pg_len division = divideBySmallInteger(pos, symbolsPerElement);
                uint_pg_len reminder = moduloBySmallInteger(pos, symbolsPerElement, division);
                return sPacker->reverseSequence(sequences[reminder], division * symbolsPerElement, length);
            };

            inline const string getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length) { 
                return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset, length); 
            };
            
            inline const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
                return orgPg + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
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
    using MultiPackedPseudoGenomeOfConstantLengthReadsType = MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
};

#endif // MULTIPACKEDPSEUDOGENOME_H_INCLUDED
