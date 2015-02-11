#ifndef SUFFIXARRAYGENERATOR_H_INCLUDED
#define SUFFIXARRAYGENERATOR_H_INCLUDED

#include "../../pseudogenome/generator/PackedPseudoGenomeGenerator.h"
#include "../DefaultSuffixArrayFactory.h"
#include "../SparseSuffixArrayFactory.h"

namespace PgSAIndex {

    class SuffixArrayGenerator {
        public:

            static SuffixArrayBase* generateDefaultSuffixArray(string readsFile, PseudoGenomeGeneratorFactory* pggf) {

                PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsFile);
                
                // FIXME: handle destruction of pgb...
                PseudoGenomeBase* pgb = pggb->generatePseudoGenomeBase();
                delete(pggb); // frees lots of memory :)...
                
                return generateSuffixArray(pgb, PGSATYPE_DEFAULT);
            }
            
            static SuffixArrayBase* generateDefaultSuffixArray(string readsFile, string pairFile, PseudoGenomeGeneratorFactory* pggf) {

                PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsFile, pairFile);
                
                // FIXME: handle destruction of pgb...
                PseudoGenomeBase* pgb = pggb->generatePseudoGenomeBase();
                delete(pggb); // frees lots of memory :)...
                
                return generateSuffixArray(pgb, PGSATYPE_DEFAULT);
            }
            
            static SuffixArrayBase* generateDefaultSuffixArray(PseudoGenomeBase* pgb) {
                return generateSuffixArray(pgb, PGSATYPE_DEFAULT);
            }
            
            static SuffixArrayBase* generateSparseSuffixArray(PseudoGenomeBase* pgb, uchar symbolsInterval) {

                PseudoGenomeGeneratorBase* pggb = new PackedPseudoGenomeGenerator(pgb, symbolsInterval);
                
                PseudoGenomeBase* ppgb = pggb->generatePseudoGenomeBase();
                
                delete(pggb);
                
                return generateSuffixArray(ppgb, PGSATYPE_SPARSE);
            }

            // WARNING: This method skips generation of duplicate filter (DONE in SA)
            static SuffixArrayBase* generateSparseSuffixArray(string readsFile, PseudoGenomeGeneratorFactory* pggf, uchar symbolsInterval) {
                PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsFile);
                
                // FIXME: handle destruction of pgb... now necessary for maintaining ReadsList and Cache
                PseudoGenomeBase* pgb = pggb->generatePseudoGenomeBase();
                delete(pggb); // frees lots of memory :
                
                return generateSparseSuffixArray(pgb, symbolsInterval);
            }

        private:
            static SuffixArrayBase* generateSuffixArray(PseudoGenomeBase* pgb, string saTypeID) {
                if (pgb->isReadLengthMin()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd()) {
                            return generateSuffixArrayTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>(pgb, saTypeID);
                        }
                        if (pgb->isPGLengthMax())
                            return generateSuffixArrayTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>(pgb, saTypeID);
                    }
                }
                if (pgb->isReadLengthStd()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return generateSuffixArrayTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>(pgb, saTypeID);
                        if (pgb->isPGLengthMax())
                            return generateSuffixArrayTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>(pgb, saTypeID);
                    }
                }

                cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
                return 0;
            }

            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static SuffixArrayBase* generateSuffixArrayTemplate(PseudoGenomeBase* pgb, string saTypeID) {               
                if (saTypeID == PGSATYPE_DEFAULT)
                    return generateDefaultSuffixArrayTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>(pgb);
                
                if (saTypeID == PGSATYPE_SPARSE) {
                    PackedPseudoGenomeBase* ppgb = static_cast<PackedPseudoGenomeBase*>(pgb);
                    if (ppgb->isPgElementMinimal()) 
                        return generateSparseSuffixArrayTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min>(ppgb);
                    if (ppgb->isPgElementStandard()) 
                        return generateSparseSuffixArrayTemplate<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std>(ppgb);

                }
                    
                
                cout << "ERROR: unknown PGSATYPE " << saTypeID;
                return 0;
            }
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static SuffixArrayBase* generateDefaultSuffixArrayTemplate(PseudoGenomeBase* pgb) {
                if (pgb->isReadLengthConstant()) {
                    typedef DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len> PseudoGenomeClass;
                    PseudoGenomeClass* pg = PseudoGenomeClass::castBase(pgb);
                    return DefaultSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len>::getSuffixArrayOfConstantLenghtReads(pg);
                }

                cout << "ERROR: unsupported variable length reads in " << PGSATYPE_DEFAULT;
                return 0;
            }
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element >
            static SuffixArrayBase* generateSparseSuffixArrayTemplate(PackedPseudoGenomeBase* ppgb) {
                const uchar skippedSymbolsCount = ppgb->getSymbolsPerElement() - 1;
                const uchar symbolsCount = ppgb->getReadsSetProperties()->symbolsCount;
                
                if (ppgb->isReadLengthConstant()) {
                    typedef PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element> PseudoGenomeClass;
                    PseudoGenomeClass* pg = PseudoGenomeClass::castBase(ppgb);
                    if(SymbolsPackingFacility<uint_ps_element_min>::isCompatibile(skippedSymbolsCount, symbolsCount))
                        return SparseSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_min>::getSuffixArrayOfConstantLenghtReads(pg);
                    if(SymbolsPackingFacility<uint_ps_element_std>::isCompatibile(skippedSymbolsCount, symbolsCount))
                        return SparseSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_std>::getSuffixArrayOfConstantLenghtReads(pg);
                }

                cout << "ERROR: unsupported variant of " << PGSATYPE_DEFAULT;
                return 0;
            }

    };

}

#endif // SUFFIXARRAYGENERATOR_H_INCLUDED
