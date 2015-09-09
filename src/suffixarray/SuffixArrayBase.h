#ifndef SUFFIXARRAYBASE_H_INCLUDED
#define SUFFIXARRAYBASE_H_INCLUDED

#include "../pseudogenome/PseudoGenomeBase.h"

namespace PgSAIndex {

    template<typename uint_pg_len>
    struct SARange {
        uint_pg_len start, end;
    };

    class SuffixArrayBase
    {
        protected:
            PseudoGenomeBase* const pgBase;
            // FIXME: should be const
            uint_max suffixArrayBytes;
            uint_max elementsCount;
            
            SuffixArrayBase(uint_max elementsCount, uint_max suffixArrayBytes, PseudoGenomeBase* pg)
            :   pgBase(pg),
                suffixArrayBytes(suffixArrayBytes),
                elementsCount(elementsCount)
            { };

        public:

            virtual ~SuffixArrayBase() {};
            
            PseudoGenomeBase* const getPseudoGenomeBase() { return pgBase; }
            
            const uint_max getSizeInBytes() {
                return this->suffixArrayBytes;             
            };
            
            inline const uint_max getElementsCount() { 
                return this->elementsCount; 
            };

            virtual string getTypeID() = 0;
            virtual void write(std::ostream& dest) = 0;

            static const string PGSA_FILE_SUFFIX;
    };

    class SuffixArrayHeader {
        private:
            string type;
            uint_read_len_max maxReadLength;
            uint_reads_cnt_max readsCount;
            uint_pg_len_max pgLength;
            uint_max saSizeInBytes;

        public:

            const static string PGSA_HEADER;

            SuffixArrayHeader(SuffixArrayBase* base) {
                this->type = base->getTypeID();
                this->maxReadLength = base->getPseudoGenomeBase()->getReadsSetProperties()->maxReadLength;
                this->readsCount = base->getPseudoGenomeBase()->getReadsSetProperties()->readsCount;
                this->pgLength = base->getPseudoGenomeBase()->getPseudoGenomeLength();
                this->saSizeInBytes = base->getSizeInBytes();
            }

            SuffixArrayHeader(std::istream& src) {

                string line;
                src >> line;
                if (line != PGSA_HEADER)
                    cout << "WARNING: wrong PGSA_HEADER";

                src >> type;

                src >> maxReadLength;
                src >> readsCount;
                src >> pgLength;
                src >> saSizeInBytes;
                src.get();
            }

            void write(std::ostream& dest) {

                dest << PGSA_HEADER << "\n";
                dest << type << "\n";

                dest << maxReadLength << "\n";
                dest << readsCount << "\n";
                dest << pgLength << "\n";
                dest << saSizeInBytes << "\n";
            }

            string getType() { return this->type; };

            bool isReadLengthMin() { return PgSAReadsSet::isReadLengthMin(maxReadLength); };
            bool isReadLengthStd() { return PgSAReadsSet::isReadLengthStd(maxReadLength); };

            bool isReadsCountStd() { return PgSAReadsSet::isReadsCountStd(readsCount); };
            bool isReadsCountMax() { return PgSAReadsSet::isReadsCountMax(readsCount); };

            bool isPGLengthStd() { return PgSAIndex::isPGLengthStd(pgLength); };
            bool isPGLengthMax() { return PgSAIndex::isPGLengthMax(pgLength); };

            bool isSALengthStd() { return PgSAIndex::isSALengthStd(pgLength); };
            bool isSALengthMax() { return PgSAIndex::isSALengthMax(pgLength); };

            uint_pg_len_max getPseudoGenomeLength() { return this->pgLength; };
    };

}

#endif // SUFFIXARRAYBASE_H_INCLUDED
