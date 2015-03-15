#ifndef DEFAULTSUFFIXARRAYLOOKUPTABLE_H_INCLUDED
#define DEFAULTSUFFIXARRAYLOOKUPTABLE_H_INCLUDED

#include "SuffixArrayLookupTableInterface.h"
#include "../../readsset/ReadsSetBase.h"
#include "../../suffixarray/SuffixArrayInterface.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

namespace PgSAIndex {

    const string SALUT_TYPE_DEFAULT = "DEFAULT_SALUT";
    
    template <  typename uint_read_len,
                typename uint_pg_len, 
                class SuffixArrayClass>
    class DefaultSuffixArrayLookupTable: public SuffixArrayLookupTableInterface<uint_read_len, uint_pg_len>
    {
        private:

            const uint_symbols_cnt symbolsCount = 0;
            const char* const symbolsList;
            const int* const symbolOrder;

            uint_pg_len* lookup;

            inline uint_max pgSuffixToLookupTableIdx(const char_pg* suffix) {
                uint_max idx = 0;
                for(uint_read_len i = 0; i < this->keyPrefixLength; i++)
                    idx = idx * symbolsCount + symbolOrder[(unsigned char) *suffix++];

                return idx;
            }

            inline uint_max pgSuffixToLookupTableIdx(const char_pg* suffix,  const uint_read_len& kmerLength, short appendSymbolOrder) {
                uint_max idx = 0;
                uint_read_len i = 0;
                while (i++ < kmerLength)
                    idx = idx * symbolsCount + symbolOrder[(unsigned char) *suffix++];
                while (i++ <= this->keyPrefixLength) // <= because previous loop additionally increased i
                    idx = idx * symbolsCount + appendSymbolOrder;
                return idx;
            }
            
        public:
            
            DefaultSuffixArrayLookupTable(uint_read_len keyPrefixLength, ReadsSetProperties* properties)
            :   SuffixArrayLookupTableInterface<uint_read_len, uint_pg_len>(keyPrefixLength),
                symbolsCount(properties->symbolsCount),
                symbolsList(properties->symbolsList),
                symbolOrder(properties->symbolOrder)
            {
            }

            ~DefaultSuffixArrayLookupTable() {
                delete[] lookup;
            };

            void read(std::istream& src)
            {
                string line;
                src >> line;

                if (line != SALUT_TYPE_DEFAULT)
                    cout << "WARNING: wrong SALUT_TYPE_DEFAULT";

                uint_max lookupTableLength;
                src >> lookupTableLength;
                src.get(); // '/n'

                if (lookupTableLength != this->getLookupTableLengthWithGuard())
                    cout << "WARNING: wrong size of lookup table.";

                lookup = (uint_pg_len*) PgSAHelpers::readArray(src, this->getLookupTableLengthWithGuard() * sizeof(uint_pg_len));
            }

            template<typename uint_reads_cnt>
            void generate(SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>* pgsa, uint_pg_len saElementsCount) {
                clock_checkpoint();
                
                lookup = new uint_pg_len[getLookupTableLengthWithGuard()];
                
                uint_max j = 0;
                for(uint_pg_len i = 0; i < saElementsCount; i++) {
                    const string suffix = pgsa->getSuffix(i, this->keyPrefixLength);
                    uint_max lutIdx = pgSuffixToLookupTableIdx(suffix.data());
                    while (j <= lutIdx)
                        lookup[j++] = i;
                }
                while(j < getLookupTableLengthWithGuard())
                    lookup[j++] = saElementsCount;
                
                cout << "SA LUT generation time " << clock_millis() << " msec!\n";
            }
            
            void write(ostream& dest) {
                dest << SALUT_TYPE_DEFAULT << "\n";

                dest << this->getLookupTableLengthWithGuard() << "\n";
                PgSAHelpers::writeArray(dest, lookup, this->getLookupTableLengthWithGuard() * sizeof(uint_pg_len));
                dest << "\n";
            };

            // returns true if the range is exact
            inline bool findSARange(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { 
            
                if (kmerLength >= this->keyPrefixLength) {
                    uint_max idx = pgSuffixToLookupTableIdx(kmer);
                    range.start = lookup[idx];
                    range.end = lookup[idx+1];
                    
                    return (kmerLength == this->keyPrefixLength);
                }
                
                range.start = lookup[pgSuffixToLookupTableIdx(kmer, kmerLength, 0)];
                range.end = lookup[pgSuffixToLookupTableIdx(kmer, kmerLength, symbolsCount-1) + 1];

                return true;                
            };
            
            // index may not be exact
            inline uint_pg_len findSAIndex(const char* kmer, const uint_read_len& kmerLength) { 
            
                if (kmerLength >= this->keyPrefixLength) {
                    uint_max idx = pgSuffixToLookupTableIdx(kmer);
                    return lookup[idx];
                }
                
                return lookup[pgSuffixToLookupTableIdx(kmer, kmerLength, 0)];
            };
            
            uint_max getLookupTableLengthWithGuard() {
                return powuint(symbolsCount, this->getKeyPrefixLength()) + 1;
            };

            uint_pg_len getRawValue(uint_max index) {
                return lookup[index];
            };

    };
}

#endif // DEFAULTSUFFIXARRAYLOOKUPTABLE_H_INCLUDED
