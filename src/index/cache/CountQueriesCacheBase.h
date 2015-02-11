/* 
 * File:   CountQueriesCacheBase.h
 * Author: Tomek
 *
 * Created on 11 luty 2014, 14:54
 */

#ifndef COUNTQUERIESCACHEBASE_H
#define	COUNTQUERIESCACHEBASE_H

#include "../../helper.h"
#include "../../pgsaconfig.h"

namespace PgSAIndex {

    struct CountQueriesResults {
        static const uint_reads_cnt_std countUnknown = UINT_MAX;
        static const uint_max occurrencesUnknown = ULLONG_MAX;
        
        uint_reads_cnt_std readsCount;
        uint_max occurrencesCount;
        uint_reads_cnt_std singleOccurrencesCount;
    };
    
    static const CountQueriesResults unknownCacheCount = { UINT_MAX, ULLONG_MAX, UINT_MAX };
    
    class CountQueriesCacheBase
    {
        protected:
            CountQueriesCacheBase()
            { };

        public:

            virtual ~CountQueriesCacheBase() {};

            virtual string getTypeID() = 0;
            virtual void write(std::ostream& dest) = 0;

            static const string CACHE_FILE_SUFFIX;
    };

    class CountQueriesCacheHeader {
        private:
            string type;

        public:

            static const string CACHE_HEADER;

            CountQueriesCacheHeader(CountQueriesCacheBase* base) {
                this->type = base->getTypeID();
            }

            CountQueriesCacheHeader(std::istream& src) {
                string line;
                src >> line;
                if (line != CACHE_HEADER)
                    cout << "WARNING: wrong CACHE_HEADER";
                src >> type;
                src.get();
            }

            void write(std::ostream& dest) {
                dest << CACHE_HEADER << "\n";
                dest << type << "\n";
            }

            string getType() { return this->type; };

    };
    
}

#endif	/* COUNTQUERIESCACHEBASE_H */

