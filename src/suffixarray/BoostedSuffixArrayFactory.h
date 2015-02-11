#ifndef BOOSTEDSUFFIXARRAYFACTORY_H_INCLUDED
#define BOOSTEDSUFFIXARRAYFACTORY_H_INCLUDED

#include "BoostedSuffixArray.h"
#include "../pseudogenome/readslist/ReadsListTypes.h"

namespace PgSAIndex {

    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                typename uint_pg_element,
                class ReadsListClass>
    class WrongBoostedSuffixArray {
        public:
            WrongBoostedSuffixArray(MultiPackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>*, std::istream& src) {
                cout << "ERROR: INCORRECT BoostedSuffixArray TEMPLATE!!!\n";
            };
    };
    
    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    using WrongBoostedSuffixArrayOfConstantLengthType = WrongBoostedSuffixArray<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, char bytesPerReadIndex >
    struct BoostedSuffixArrayOfConstantLengthTypeTemplate {
        typedef WrongBoostedSuffixArrayOfConstantLengthType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element> Type;
    };

    template<typename uint_pg_len, typename uint_pg_element>
    struct BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 3> {
        typedef BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 4, 3, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element>
    struct BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 4> {
        typedef BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 5, 4, 0xFFFFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element>
    struct BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 3> {
        typedef BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 5, 3, 0x00FFFFFF> Type;
    };
    
    template<typename uint_pg_len, typename uint_pg_element>
    struct BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 4> {
        typedef BoostedSuffixArrayOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len, uint_pg_element, 6, 4, 0xFFFFFFFF> Type;
    };
 
    template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    class BoostedSuffixArrayFactory {
        public:

            static SuffixArrayBase* getSuffixArrayOfConstantLenghtReads(MultiPackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element>* mppg, std::istream& src) {
                if (bytesPerValue(mppg->readsCount()) <= BytesPerReadIndex<uint_reads_cnt>::minimum)
                    return new typename BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type(mppg, src);
                if (bytesPerValue(mppg->readsCount()) == BytesPerReadIndex<uint_reads_cnt>::standard)
                    return new typename BoostedSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, BytesPerReadIndex<uint_reads_cnt>::standard>::Type(mppg, src);
                cout << "WARNING: SA factory fail!";
                return 0;
            }
    };

}

#endif // BOOSTEDSUFFIXARRAYFACTORY_H_INCLUDED
