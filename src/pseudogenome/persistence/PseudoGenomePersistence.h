#ifndef PSEUDOGENOMEPERSISTENCE_H
#define PSEUDOGENOMEPERSISTENCE_H

#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../../helper.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    class PseudoGenomePersistence
    {
        private:
            PseudoGenomePersistence() {};

        public:
            virtual ~PseudoGenomePersistence() {};

            static void writePseudoGenome(PseudoGenomeBase* pgb, std::ostream& dest) {
                // header
                PseudoGenomeHeader pgh(pgb);                
                pgh.write(dest);
                if (pgh.getType() == PGTYPE_PACKED) {
                    PackedPseudoGenomeHeaderExtension ppghe(static_cast<PackedPseudoGenomeBase*>(pgb));
                    ppghe.write(dest);
                }

                // default pseudogenome
                pgb->write(dest);

            }

            static void writePseudoGenome(PseudoGenomeBase* pgb, string pseudoGenomePrefix) {
                std::ofstream dest(pseudoGenomePrefix + PseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX, std::ios::out | std::ios::binary);
                writePseudoGenome(pgb, dest);
                dest.close();
            }

            static bool isValidPseudoGenome(string file) {
                unsigned char suffixLength = PseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX.length();
                return ((file.length() > suffixLength) && 
                    (file.substr(file.length() - suffixLength).compare(PseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX) == 0));
            }
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static PseudoGenomeBase* readPseudoGenomeTemplate(std::istream& src, PseudoGenomeHeader& pgh);

            static PseudoGenomeBase* readPseudoGenome(std::istream& src){
                PseudoGenomeHeader pgh(src);

                if (pgh.isReadLengthMin()) {
                    if (pgh.isReadsCountStd()) {
                        if (pgh.isPGLengthStd())
                            return readPseudoGenomeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>(src, pgh);
                        if (pgh.isPGLengthMax())
                            return readPseudoGenomeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>(src, pgh);
                    }
                }
                if (pgh.isReadLengthStd()) {
                    if (pgh.isReadsCountStd()) {
                        if (pgh.isPGLengthStd())
                            return readPseudoGenomeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>(src, pgh);
                        if (pgh.isPGLengthMax())
                            return readPseudoGenomeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>(src, pgh);
                    }
                }

                cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
                return 0;
            }

            static PseudoGenomeBase* readPseudoGenome(string pseudoGenomeFile){
                std::ifstream src(pseudoGenomeFile, std::ios::in | std::ios::binary);

                return readPseudoGenome(src);
            }

    };

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
    PseudoGenomeBase* PseudoGenomePersistence::readPseudoGenomeTemplate(std::istream& src, PseudoGenomeHeader& pgh) {
        if (pgh.getType() == PGTYPE_DEFAULT) {
            if (pgh.isReadLengthConstant())
                return new DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>(pgh.getPseudoGenomeLength(), src);
        }
        
        if (pgh.getType() == PGTYPE_PACKED) {
            PackedPseudoGenomeHeaderExtension ppghe(src);
            if (pgh.isReadLengthConstant()) {
                if (ppghe.isPgElementMinimal())
                    return new PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min>(pgh.getPseudoGenomeLength(), src);
                if (ppghe.isPgElementStandard())
                    return new PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std>(pgh.getPseudoGenomeLength(), src);
            }
        }
        
        cout << "ERROR: unsupported PGTYPE " << pgh.getType() << " with " << (pgh.isReadLengthConstant()?"constant":"variable") << " read length\n";
        return 0;
    }
}

#endif // PSEUDOGENOMEPERSISTENCE_H
