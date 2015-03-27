#ifndef PGSAINDEXFACTORY_H_INCLUDED
#define PGSAINDEXFACTORY_H_INCLUDED

#include "DefaultRevCompPgSAIndex.h"
#include "../suffixarray/SparseSuffixArrayBase.h"
#include "../suffixarray/SuffixArrayInterface.h"
#include "../suffixarray/persistence/SuffixArrayPersistence.h"
#include "cache/DefaultCountQueriesCache.h"
#include "cache/persistence/CountQueriesCachePersistence.h"

//TODO: Why SuffixArrayPersistance does not need to be included???

namespace PgSAIndex {

    class PgSAIndexFactory {
        private:

            template < typename api_uint_reads_cnt, typename api_uint_read_len >
            static PgSAIndexInterface<api_uint_reads_cnt, api_uint_read_len>* getPgSAIndexAPITemplate(string pgsaFile, string cacheFile, bool reverseComplements) {
                SuffixArrayBase* sab = SuffixArrayPersistence::readPgSA(pgsaFile);
                
                DefaultCountQueriesCacheStandard* dcqcs = 0;
                if (cacheFile.compare("")) 
                    try {
                        dcqcs = static_cast<DefaultCountQueriesCacheStandard*>(DefaultCountQueriesCacheStandard::castBase(
                                CountQueriesCachePersistence::readCache(cacheFile)
                                ));
                    } catch (int e) {
                    }
                
                PseudoGenomeBase* pgb = sab->getPseudoGenomeBase();

                 if (pgb->isReadLengthMin()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return getPgSAIndexBasicTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>(sab, dcqcs, pgsaFile, reverseComplements);
                        if (pgb->isPGLengthMax())
                            return getPgSAIndexBasicTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>(sab, dcqcs, pgsaFile, reverseComplements);
                    }
                }
                if (pgb->isReadLengthStd()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return getPgSAIndexBasicTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>(sab, dcqcs, pgsaFile, reverseComplements);
                        if (pgb->isPGLengthMax())
                            return getPgSAIndexBasicTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>(sab, dcqcs, pgsaFile, reverseComplements);
                    }
                }

                cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
                return 0;

            }
            
            template< typename api_uint_reads_cnt, typename api_uint_read_len, typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static PgSAIndexInterface<api_uint_reads_cnt, api_uint_read_len>* getPgSAIndexBasicTemplate(SuffixArrayBase* sab, DefaultCountQueriesCacheStandard* dcqcs, string pgsaPrefix, bool reverseComplements) {
                ReadsSetProperties* rsp = sab->getPseudoGenomeBase()->getReadsSetProperties();
                
                if (sab->getTypeID() == PGSATYPE_DEFAULT) {
                    if(rsp->constantReadLength) {
                        if (bytesPerValue(rsp->readsCount) <= BytesPerReadIndex<uint_reads_cnt>::minimum) {
                            typedef typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type SuffixArrayClass;
                            return getPgSAIndexConstructorTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>(
                                    sab, dcqcs, pgsaPrefix, reverseComplements);
                        }
                        if (bytesPerValue(rsp->readsCount) == BytesPerReadIndex<uint_reads_cnt>::standard) {
                            typedef typename DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, BytesPerReadIndex<uint_reads_cnt>::standard>::Type SuffixArrayClass;
                            return getPgSAIndexConstructorTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>(
                                    sab, dcqcs, pgsaPrefix, reverseComplements);
                        }
                    }
                }
                
                if (sab->getTypeID() == PGSATYPE_SPARSE) {
                   PackedPseudoGenomeBase* ppgb = static_cast<PackedPseudoGenomeBase*>(sab->getPseudoGenomeBase());
                   if (ppgb->isPgElementMinimal())
                        return getPgSAIndexSparseTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min>(sab, dcqcs, ppgb, pgsaPrefix, reverseComplements);
                   if (ppgb->isPgElementStandard())
                        return getPgSAIndexSparseTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std>(sab, dcqcs, ppgb, pgsaPrefix, reverseComplements);
                }
                
                cout << "ERROR: unsupported PgSA " << sab->getTypeID();
                return 0;                
            }

            template< typename api_uint_reads_cnt, typename api_uint_read_len, typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element >
            static PgSAIndexInterface<api_uint_reads_cnt, api_uint_read_len>* getPgSAIndexSparseTemplate(SuffixArrayBase* sab, DefaultCountQueriesCacheStandard* dcqcs, PackedPseudoGenomeBase* ppgb, string pgsaPrefix, bool reverseComplements) {
                SparseSuffixArrayBase* ssab = static_cast<SparseSuffixArrayBase*>(sab);
                if (ssab->isSAElementMinimal())
                    return getPgSAIndexSparseFullTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_min>(sab, dcqcs, ppgb, pgsaPrefix, reverseComplements);
                if (ssab->isSAElementStandard())
                    return getPgSAIndexSparseFullTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_std>(sab, dcqcs, ppgb, pgsaPrefix, reverseComplements);
               
                cout << "ERROR: unsupported PgSA " << sab->getTypeID();
                return 0;             
            }
                        
            template< typename api_uint_reads_cnt, typename api_uint_read_len, typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, typename uint_skipped_element>
            static PgSAIndexInterface<api_uint_reads_cnt, api_uint_read_len>* getPgSAIndexSparseFullTemplate(SuffixArrayBase* sab, DefaultCountQueriesCacheStandard* dcqcs, PackedPseudoGenomeBase* ppgb, string pgsaPrefix, bool reverseComplements) {
                if (ppgb->isReadLengthConstant()) {
                    if (bytesPerValue(ppgb->getReadsSetProperties()->readsCount) <= BytesPerReadIndex<uint_reads_cnt>::minimum) {
                        typedef typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::minimum>::Type SuffixArrayClass;
                        return getPgSAIndexConstructorTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>(
                                        sab, dcqcs, pgsaPrefix, reverseComplements);
                    }
                    if (bytesPerValue(ppgb->getReadsSetProperties()->readsCount) <= BytesPerReadIndex<uint_reads_cnt>::standard) {
                        typedef typename SparseSuffixArrayOfConstantLengthTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_skipped_element, BytesPerReadIndex<uint_reads_cnt>::standard>::Type SuffixArrayClass;
                        return getPgSAIndexConstructorTemplate<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>(
                                        sab, dcqcs, pgsaPrefix, reverseComplements);
                    }
                }    
                cout << "ERROR: unsupported PgSA " << sab->getTypeID();
                return 0;             
            }

            template< typename api_uint_reads_cnt, typename api_uint_read_len, typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class SuffixArrayClass>
            static PgSAIndexInterface<api_uint_reads_cnt, api_uint_read_len>* getPgSAIndexConstructorTemplate(SuffixArrayBase* sab, DefaultCountQueriesCacheStandard* dcqcs, string pgsaPrefix, bool reverseComplements) {
                typedef DefaultPgSAIndex<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass> DefaultPgSAIndexClass;
                DefaultPgSAIndexClass* idx = new DefaultPgSAIndexClass(*SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>::castBase(sab), dcqcs, pgsaPrefix);
                if (reverseComplements) {
                    typedef DefaultRevCompPgSAIndex<api_uint_reads_cnt, api_uint_read_len, uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass> DefaultRevCompPgSAIndexClass;
                    return new DefaultRevCompPgSAIndexClass(idx);
                }else 
                    return idx;
            }
            
        public:
            static PgSAIndexStandard* getPgSAIndexStandard(string pgsaFile, string cacheFile, bool reverseComplements);
            
            static PgSAIndexStandard* getPgSAIndexStandard(string pgsaFile, bool boosted, bool reverseComplements) {
                return getPgSAIndexAPITemplate<uint_reads_cnt_std, unsigned int >(pgsaFile, "", reverseComplements);
            }
    };
    
}

#endif // PGSAINDEXFACTORY_H_INCLUDED
