#ifndef SUFFIXARRAYPERSISTENCE_H
#define SUFFIXARRAYPERSISTENCE_H

#include "../../suffixarray/DefaultSuffixArrayFactory.h"
#include "../../suffixarray/SparseSuffixArrayFactory.h"
#include "../../pseudogenome/persistence/PseudoGenomePersistence.h"
#include "../../helper.h"

namespace PgSAIndex {

    class SuffixArrayPersistence
    {
        private:
            SuffixArrayPersistence() {};

        public:
            virtual ~SuffixArrayPersistence() {};

            const static string DEFAULT_PGSATYPE;

            static void writePgSA(SuffixArrayBase* pgsab, string pgsaPrefix) {

                std::ofstream dest(pgsaPrefix + SuffixArrayBase::PGSA_FILE_SUFFIX, std::ios::out | std::ios::binary);
                
                SuffixArrayHeader pgsah(pgsab);
                pgsah.write(dest);            

                PseudoGenomePersistence::writePseudoGenome(pgsab->getPseudoGenomeBase(), dest);

                if (pgsah.getType() == PGSATYPE_SPARSE) {
                    SparseSuffixArrayHeaderExtension ssahe(static_cast<SparseSuffixArrayBase*>(pgsab));
                    ssahe.write(dest);
                }
                
                pgsab->write(dest);
                dest.close();
            };

            static bool isValidSuffixArray(string file) {
                unsigned char len = SuffixArrayBase::PGSA_FILE_SUFFIX.length();
                return ((file.length() > len) && 
                    (file.substr(file.length() - len).compare(SuffixArrayBase::PGSA_FILE_SUFFIX)) == 0);
            }
            
            static PseudoGenomeBase* readPgOnly(string pgsaFile) {
                std::ifstream src(pgsaFile, std::ios::in | std::ios::binary);

                SuffixArrayHeader pgsah(src);
                
                return PseudoGenomePersistence::readPseudoGenome(src);
            }
            
            static SuffixArrayBase* readPgSA(string pgsaFile) {
                std::ifstream src(pgsaFile, std::ios::in | std::ios::binary);

                SuffixArrayHeader pgsah(src);
                
                PseudoGenomeBase* pgb = PseudoGenomePersistence::readPseudoGenome(src);
                 if (pgb->isReadLengthMin()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return readPgSA<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>(src, pgb, pgsah);
                        if (pgb->isPGLengthMax())
                            return readPgSA<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>(src, pgb, pgsah);
                    }
                }
                if (pgb->isReadLengthStd()) {
                    if (pgb->isReadsCountStd()) {
                        if (pgb->isPGLengthStd())
                            return readPgSA<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>(src, pgb, pgsah);
                        if (pgb->isPGLengthMax())
                            return readPgSA<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>(src, pgb, pgsah);
                    }
                }

                cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
                return 0;
            };

            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static SuffixArrayBase* readPgSA(std::istream& src, PseudoGenomeBase* pgb, SuffixArrayHeader& pgsah) {

                if (pgsah.getType() == PGSATYPE_DEFAULT)
                    if(pgb->isReadLengthConstant())
                        return DefaultSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len>::getSuffixArrayOfConstantLenghtReads(DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>::castBase(pgb), src);

                 
                if (pgsah.getType() == PGSATYPE_SPARSE) {
                    PackedPseudoGenomeBase* ppgb = static_cast<PackedPseudoGenomeBase*>(pgb);
                    if (ppgb->isPgElementMinimal()) 
                        return readSparsePgSA<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min>(src, ppgb);
                    if (ppgb->isPgElementStandard()) 
                        return readSparsePgSA<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std>(src, ppgb);
                }
                 
                cout << "ERROR: unknown PGSATYPE " << pgsah.getType();
                return 0;
            };
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element >
            static SuffixArrayBase* readSparsePgSA(std::istream& src, PackedPseudoGenomeBase* ppgb) {
                SparseSuffixArrayHeaderExtension ssahe(src);
                
                if (ppgb->isReadLengthConstant()) {
                    typedef PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element> PseudoGenomeClass;
                    PseudoGenomeClass* pg = PseudoGenomeClass::castBase(ppgb);
                    if(ssahe.isSAElementMinimal())
                        return SparseSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_min>::getSuffixArrayOfConstantLenghtReads(pg, src);
                    if(ssahe.isSAElementStandard())
                        return SparseSuffixArrayFactory<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, uint_ps_element_std>::getSuffixArrayOfConstantLenghtReads(pg, src);
                }
                 
                cout << "ERROR: unsupported variant of sparse SA";
                return 0;
            };
    };
}

#endif // SUFFIXARRAYPERSISTENCE_H
