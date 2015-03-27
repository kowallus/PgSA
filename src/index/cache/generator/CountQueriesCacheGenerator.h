/* 
 * File:   CountQueriesCacheGenerator.h
 * Author: coach
 *
 * Created on 26 marca 2015, 22:53
 */

#ifndef COUNTQUERIESCACHEGENERATOR_H
#define	COUNTQUERIESCACHEGENERATOR_H

#include "../DefaultCountQueriesCache.h"
#include "../../../pseudogenome/DefaultPseudoGenome.h"



namespace PgSAIndex {

    class CountQueriesCacheGenerator {
        public:

            static CountQueriesCacheBase* generateCountQueriesCache(PseudoGenomeBase* pgb) {
                return generateCountQueriesCache(pgb, PGTYPE_DEFAULT);
            }


        private:
            static CountQueriesCacheBase* generateCountQueriesCache(PseudoGenomeBase* pgb, string pgTypeID) {
                if (pgb->isReadLengthMin()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd()) {
                            return generateCountQueriesCacheTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>(pgb, pgTypeID);
                        }
                        if (pgb->isPGLengthMax())
                            return generateCountQueriesCacheTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>(pgb, pgTypeID);
                    }
                }
                if (pgb->isReadLengthStd()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return generateCountQueriesCacheTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>(pgb, pgTypeID);
                        if (pgb->isPGLengthMax())
                            return generateCountQueriesCacheTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>(pgb, pgTypeID);
                    }
                }

                cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
                return 0;
            }

            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static CountQueriesCacheBase* generateCountQueriesCacheTemplate(PseudoGenomeBase* pgb, string pgTypeID) {               
                if (pgTypeID == PGTYPE_DEFAULT)
                    return generateCountQueriesCacheTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>(pgb);
                
                cout << "ERROR: unknown PGSATYPE " << pgTypeID;
                return 0;
            }
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static CountQueriesCacheBase* generateCountQueriesCacheTemplate(PseudoGenomeBase* pgb) {
                if (pgb->isReadLengthConstant()) {
                    typedef DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len> PseudoGenomeClass;
                    PseudoGenomeClass* pg = PseudoGenomeClass::castBase(pgb);
                    return new GeneratedCountQueriesCache(pg, pg->getReadsSetProperties());
                }

                cout << "ERROR: unsupported variable length reads in " << PGTYPE_DEFAULT;
                return 0;
            }
    };

}



#endif	/* COUNTQUERIESCACHEGENERATOR_H */

