#ifndef SPARSESUFFIXARRAYFACTORY_H_INCLUDED
#define SPARSESUFFIXARRAYFACTORY_H_INCLUDED

#include "SparseSuffixArray.h"
#include "../pseudogenome/readslist/ReadsListTypes.h"

namespace PgSAIndex {

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len, typename uint_pg_element,
                class ReadsListClass>
    class WrongSparseSuffixArray {
        public:
            WrongSparseSuffixArray(PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* pseudoGenome) {
                cout << "ERROR: INCORRECT SparseSuffixArray TEMPLATE!!!\n";
            };
            WrongSparseSuffixArray(PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* pseudoGenome, std::istream& src) {
                cout << "ERROR: INCORRECT SparseSuffixArray TEMPLATE!!!\n";
            };
    };
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    using WrongSparseSuffixArrayOfConstantLengthType = WrongSparseSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element, char bytesPerReadIndex>
    struct SparseSuffixArrayOfConstantLengthTypeTemplate {
        typedef WrongSparseSuffixArrayOfConstantLengthType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element> Type;
    };

    template<typename uint_pg_len, typename uint_pg_element >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min,  3> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 5, 3, 4, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 4> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 6, 4, 5, 0xFFFFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element  >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 3> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 6, 3, 5, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element  >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 4> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_min, 7, 4, 6, 0xFFFFFFFF> Type;
    };

    template<typename uint_pg_len, typename uint_pg_element >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std,  3> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 6, 3, 4, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 4> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 7, 4, 5, 0xFFFFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element  >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 3> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 7, 3, 5, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element  >
    struct SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 4> {
        typedef SparseSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, uint_ps_element_std, 8, 4, 6, 0xFFFFFFFF> Type;
    };
    
    template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element>
    class SparseSuffixArrayFactory {       
        public:

            static SuffixArrayBase* getSuffixArrayOfConstantLenghtReads(PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element>* pg, int fixed_min_k = 1) {
                if (bytesPerValue(pg->readsCount()) <= BytesPerReadIndex<uint_reads_cnt>::minimum) 
                    return new typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type(pg, fixed_min_k);
                if (bytesPerValue(pg->readsCount()) == BytesPerReadIndex<uint_reads_cnt>::standard)
                    return new typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::standard>::Type(pg, fixed_min_k);
                cout << "WARNING: SA factory fail!";
                return 0;
            }

            static SuffixArrayBase* getSuffixArrayOfConstantLenghtReads(PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element>* pg, std::istream& src) {
                if (bytesPerValue(pg->readsCount()) <= BytesPerReadIndex<uint_reads_cnt>::minimum)
                    return new typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type(pg, src);
                if (bytesPerValue(pg->readsCount()) == BytesPerReadIndex<uint_reads_cnt>::standard)
                    return new typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::standard>::Type(pg, src);
                cout << "WARNING: SA factory fail!";
                return 0;
            }
    };

}

#endif // SPARSESUFFIXARRAYFACTORY_H_INCLUDED
