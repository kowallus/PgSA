#ifndef PSEUDOGENOMEINTERFACE_H_INCLUDED
#define PSEUDOGENOMEINTERFACE_H_INCLUDED

#include "../readsset/ReadsSetInterface.h"
#include "../readsset/ReadsSetBase.h"
#include "readslist/ReadsListInterface.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace PgSAIndex {

    //FIXME: Probably unnecessary interface... maybe apply CRTP?
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
    class PseudoGenomeInterface: public ReadsSetInterface<uint_read_len, uint_reads_cnt>
    {
        public:
            virtual ~PseudoGenomeInterface() {};
    };

    typedef PseudoGenomeInterface< uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std > MinimalPseudoGenomeInterface;
    typedef PseudoGenomeInterface< uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std > StandardPseudoGenomeInterface;
    typedef PseudoGenomeInterface< uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max > ExtendedPseudoGenomeInterface;

}

#endif // PSEUDOGENOMEINTERFACE_H_INCLUDED
