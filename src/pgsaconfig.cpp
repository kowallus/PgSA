#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "suffixarray/persistence/SuffixArrayPersistence.h"
#include "index/cache/persistence/CountQueriesCachePersistence.h"

namespace PgSAIndex {

    const string PseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX = ".pgen";

    const string PseudoGenomeHeader::PSEUDOGENOME_HEADER = "PGEN";
    
    const string SuffixArrayBase::PGSA_FILE_SUFFIX = ".pgsa";

    const string SuffixArrayHeader::PGSA_HEADER = "PgSA";

    const string CountQueriesCacheBase::CACHE_FILE_SUFFIX = ".pgc";
    
    const string CountQueriesCacheHeader::CACHE_HEADER = "PgCache";
    
}

