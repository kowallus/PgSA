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
    class MultiPackedPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>>,
            public PackedPseudoGenomeBase
    {
        private:
            uint_pg_element** sequences;
            const char_pg* orgPg;
            
            SymbolsPackingFacility<uint_pg_element>* sPacker;
            
            PseudoGenomeBase* srcPGB; //manages readsList;
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

            CountQueriesCacheBase* countQueriesCache = 0;
            
            const uint_pg_element* getRawSuffix(const uint_pg_len pos);
            
        public:

            MultiPackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement);

            ~MultiPackedPseudoGenome();

            static MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* castBase(PseudoGenomeBase* base);

            string getTypeID();

            void write(std::ostream& dest);
            
            CountQueriesCacheBase* getCountQueriesCacheBase();
            
            uint_pg_len getElementsCountWithGuard();
            
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList();

//            inline uint_pg_len getLength() { return this->length; };

            char getSymbol(uint_pg_len pos);

            // calculated position should be divisible by symbolsPerElement
            const uint_pg_element* getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos);
            
            uint_read_len maxReadLength();

            uint_reads_cnt readsCount();

            bool isReadLengthConstant();

            const string getSuffix(const uint_pg_len pos, const uint_pg_len length);

            const string getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length);
            
            const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos);
            
            const string getRead(uint_reads_cnt originalIdx);
            
            uint_read_len readLength(uint_reads_cnt originalIdx);

            const char_pg getSymbolImpl(const uint_pg_len posIdx);
            const uint_pg_len getLengthImpl();
            
            SymbolsPackingFacility<uint_pg_element>* getSymbolsPacker();
            
            uint_read_len maxReadLengthVirtual();
            uint_reads_cnt readsCountVirtual();
            bool isReadLengthConstantVirtual();
            const string getReadVirtual(uint_reads_cnt i);
            uint_read_len readLengthVirtual(uint_reads_cnt i);
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    using MultiPackedPseudoGenomeOfConstantLengthReadsType = MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
};

#endif // MULTIPACKEDPSEUDOGENOME_H_INCLUDED
