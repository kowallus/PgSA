#include "PgSAIndexFactory.h"

namespace PgSAIndex {
    
    PgSAIndexStandard* PgSAIndexFactory::getPgSAIndexStandard(string pgsaFile, string cacheFile, bool reverseComplements) {
        return getPgSAIndexAPITemplate<uint_reads_cnt_std, unsigned int >(pgsaFile, cacheFile, reverseComplements);
    }
    
}
