#ifndef DEFAULTSUFFIXARRAYFACTORY_H_INCLUDED
#define DEFAULTSUFFIXARRAYFACTORY_H_INCLUDED

#include "DefaultSuffixArray.h"
#include "../pseudogenome/readslist/ReadsListTypes.h"

namespace PgSAIndex {

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                class ReadsListClass>
    class WrongDefaultSuffixArray {
        public:
            WrongDefaultSuffixArray(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* pseudoGenome) {
                cout << "ERROR: INCORRECT DefaultSuffixArray TEMPLATE!!!\n";
            };
            WrongDefaultSuffixArray(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* pseudoGenome, std::istream& src) {
                cout << "ERROR: INCORRECT DefaultSuffixArray TEMPLATE!!!\n";
            };
    };
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using WrongDefaultSuffixArrayOfConstantLengthType = WrongDefaultSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, char bytesPerReadIndex >
    struct DefaultSuffixArrayOfConstantLengthTypeTemplate {
        typedef WrongDefaultSuffixArrayOfConstantLengthType<uint_read_len, uint_reads_cnt, uint_pg_len> Type;
    };

    template<typename uint_pg_len >
    struct DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, 3> {
        typedef DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, 4, 3, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len >
    struct DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, 4> {
        typedef DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, 5, 4, 0xFFFFFFFF> Type;
    };
    
    template<typename uint_pg_len >
    struct DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, 3> {
        typedef DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, 5, 3, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len >
    struct DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, 4> {
        typedef DefaultSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, 6, 4, 0xFFFFFFFF> Type;
    };
 
//        typedef DefaultSuffixArray<uint_read_len_min, uint_reads_cnt_max, uint_pg_len, 6, 5, 0x000000FFFFFFFFFF> ExtendedMinDefaultSuffixArray;
//        typedef DefaultSuffixArray<uint_read_len_std, uint_reads_cnt_max, uint_pg_len, 7, 5, 0x000000FFFFFFFFFF> ExtendedStdDefaultSuffixArray;

    template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    class DefaultSuffixArrayFactory {
        public:

            static SuffixArrayBase* getSuffixArrayOfConstantLenghtReads(DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>* pg, int fixed_min_k = 1) {   
                if (bytesPerValue(pg->readsCount()) <= BytesPerReadIndex<uint_reads_cnt>::minimum)
                    return new typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type(pg, fixed_min_k);
                if (bytesPerValue(pg->readsCount()) == BytesPerReadIndex<uint_reads_cnt>::standard)
                    return new typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::standard>::Type(pg, fixed_min_k);
                cout << "WARNING: SA factory fail!";
                return 0;
            }

            static SuffixArrayBase* getSuffixArrayOfConstantLenghtReads(DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>* pg, std::istream& src) {
                if (bytesPerValue(pg->readsCount()) <= BytesPerReadIndex<uint_reads_cnt>::minimum)
                    return new typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type(pg, src);
                if (bytesPerValue(pg->readsCount()) == BytesPerReadIndex<uint_reads_cnt>::standard)
                    return new typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::standard>::Type(pg, src);
                cout << "WARNING: SA factory fail!";
                return 0;
            }
    };

}

#endif // DEFAULTSUFFIXARRAYFACTORY_H_INCLUDED
