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

            DefaultPseudoGenome(uint_pg_len pgLength, ReadsSetProperties* properties);
            DefaultPseudoGenome(uint_pg_len pgLength, std::istream& src);
            ~DefaultPseudoGenome();
            
            static DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* castBase(PseudoGenomeBase* base);

            string getTypeID();;

            void write(std::ostream& dest);
            
            CountQueriesCacheBase* getCountQueriesCacheBase();
            
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList();
//            uint_pg_len getLength() { return this->length; };

            uint_pg_len getLengthWithGuard();
            
            char getSymbol(uint_pg_len pos);

            const char_pg* getSuffix(const uint_pg_len pos);
            
            const char_pg* getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset);

            const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos);
            
            uint_read_len maxReadLength();

            uint_reads_cnt readsCount();

            bool isReadLengthConstant();

            const string getRead(const uint_reads_cnt originalIdx);
            
            uint_read_len readLength(const uint_reads_cnt originalIdx);

            uint_read_len maxReadLengthVirtual();
            uint_reads_cnt readsCountVirtual();
            bool isReadLengthConstantVirtual();
            const string getReadVirtual(uint_reads_cnt i);
            uint_read_len readLengthVirtual(uint_reads_cnt i);

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

            GeneratedPseudoGenome(uint_pg_len sequenceLength, ReadsSetProperties* properties);

            ~GeneratedPseudoGenome();

            void append(const string& read, uint_read_len length, uint_read_len overlap, uint_reads_cnt orgIdx);

            void validate();
            
            void setCountQueriesCache(CountQueriesCacheBase* cqcb);
            
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using GeneratedPseudoGenomeOfConstantLengthReadsType = GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

};

#endif // DEFAULTPSEUDOGENOME_H_INCLUDED
