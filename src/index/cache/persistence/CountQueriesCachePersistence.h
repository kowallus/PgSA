/* 
 * File:   CountQueriesPersistence.h
 * Author: Tomek
 *
 * Created on 11 luty 2014, 15:14
 */

#ifndef COUNTQUERIESCACHEPERSISTENCE_H
#define	COUNTQUERIESCACHEPERSISTENCE_H

#include "../../../helper.h"
#include "../DefaultCountQueriesCache.h"

namespace PgSAIndex {

    class CountQueriesCachePersistence
    {
        private:
            CountQueriesCachePersistence() {};

        public:
            virtual ~CountQueriesCachePersistence() {};

            static void writeCache(CountQueriesCacheBase* cqcb, string cachePrefix) {

                std::ofstream dest(cachePrefix + CountQueriesCacheBase::CACHE_FILE_SUFFIX, std::ios::out | std::ios::binary);

                // header
                CountQueriesCacheHeader cqch(cqcb);
                cqch.write(dest);

                cqcb->write(dest);
                dest.close();
                
            };

            static CountQueriesCacheBase* readCache(string cacheFile) {
                std::ifstream src(cacheFile, std::ios::in | std::ios::binary);

                CountQueriesCacheHeader cqch(src);

                 if (cqch.getType() == CACHETYPE_DEFAULT)
                     return new DefaultCountQueriesCacheStandard(src);
                 

                cout << "ERROR: unknown CACHETYPE " << cqch.getType();
                return 0;
            };
    };

}

#endif	/* COUNTQUERIESCACHEPERSISTENCE_H */

