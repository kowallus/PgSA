#include "PgSAIndexFactory.h"

namespace PgSAIndex {
    
    PgSAIndexStandard* PgSAIndexFactory::getPgSAIndexStandard(string pgsaFile, string cacheFile, bool boosted, bool boosized, bool reverseComplements) {
        return getPgSAIndexAPITemplate<uint_reads_cnt_std, unsigned int >(pgsaFile, cacheFile, boosted, boosized, reverseComplements);
    }
    
}
