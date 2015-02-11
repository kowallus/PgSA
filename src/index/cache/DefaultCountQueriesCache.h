/* 
 * File:   DefaultCountQueriesCache.h
 * Author: Tomek
 *
 * Created on 11 luty 2014, 19:15
 */

#ifndef DEFAULTCOUNTQUERIESCACHE_H
#define	DEFAULTCOUNTQUERIESCACHE_H

#include "CountQueriesCacheInterface.h"
#include "../../pseudogenome/DefaultPseudoGenome.h"

using namespace PgSAHelpers;

namespace PgSAIndex {

    const string CACHETYPE_DEFAULT = "DEFAULT_CACHE";
    
    template<typename api_uint_reads_cnt>
    class DefaultCountQueriesCache: public CountQueriesCacheInterface<api_uint_reads_cnt, DefaultCountQueriesCache<api_uint_reads_cnt>>    
    {
        protected:
        
            uint_read_len_std cqCacheMaxKLength = 13;
            
            uint_read_len_std shortCacheKLength = 11;
            uint_read_len_std byteCacheKLength = 13;            
            
            static const unsigned short exceedsUShortCache = 65535;
            static const unsigned char exceedsUCharCache = 255;
            
            CountQueriesResults** cqCache = 0;
            CountQueriesResults* flatCqCache = 0;

            unsigned short** cqShortCache = 0;
            unsigned short* flatCqShortCache = 0;

            unsigned char** cqByteCache = 0;
            unsigned char* flatCqByteCache = 0;

            uint_symbols_cnt cqSymbolsCount = 0;
            char symbolsList[UCHAR_MAX];
            int symbolOrder[UCHAR_MAX];
            
            void generateSymbolOrder() {
                for (ushort i = 0; i < UCHAR_MAX; i++)
                    symbolOrder[i] = -1;
                    
                for (uint_symbols_cnt i = 0; i < cqSymbolsCount; i++)
                    symbolOrder[(unsigned char) symbolsList[(unsigned char) i]] = i;
            }
            
            uint_max sizeOfCqCacheFor(uint_read_len_std kLength) {
                return powuint(cqSymbolsCount, kLength);
            }

            uint_max sizeOfFlatCqCache(uint_read_len_std minKLength, uint_read_len_std maxKLength) {
                uint_max size = 0;
                for(uint_read_len_std i = minKLength; i <= maxKLength; i++)
                    size += sizeOfCqCacheFor(i);

                return size;
            }

            int sizeOfFlatCqCache() { return sizeOfFlatCqCache(1, shortCacheKLength - 1); }

            int sizeOfFlatCqShortCache() { return sizeOfFlatCqCache(shortCacheKLength, byteCacheKLength - 1); }

            int sizeOfFlatCqByteCache() { return sizeOfFlatCqCache(byteCacheKLength, cqCacheMaxKLength); }

            void buildArrays() {
                cqCache = new CountQueriesResults*[shortCacheKLength];
                CountQueriesResults* cqPtr = flatCqCache;
                for(uint_read_len_std i = 1; i < shortCacheKLength; i++) {
                    cqCache[i] = cqPtr;
                    cqPtr += sizeOfCqCacheFor(i);
                }

                cqShortCache = new unsigned short*[byteCacheKLength];
                unsigned short* cqShortPtr = flatCqShortCache;
                for(uint_read_len_std i = shortCacheKLength; i < byteCacheKLength; i++) {
                    cqShortCache[i] = cqShortPtr;
                    cqShortPtr += sizeOfCqCacheFor(i);
                }

                cqByteCache = new unsigned char*[(cqCacheMaxKLength + 1)];
                unsigned char* cqBytePtr = flatCqByteCache;
                for(uint_read_len_std i = byteCacheKLength; i <= cqCacheMaxKLength; i++) {
                    cqByteCache[i] = cqBytePtr;
                    cqBytePtr += sizeOfCqCacheFor(i);
                }
            }
            
            inline CountQueriesResults convert(unsigned char res) {
                CountQueriesResults cqr;
                if (res == exceedsUCharCache)
                    return unknownCacheCount;
                cqr.readsCount = res;
                cqr.occurrencesCount = res;
                cqr.singleOccurrencesCount = res;
                return cqr;
            }

            inline CountQueriesResults convert(unsigned short res) {
                CountQueriesResults cqr;
                if (res == exceedsUShortCache)
                    return unknownCacheCount;
                cqr.readsCount = res;
                cqr.occurrencesCount = res;
                cqr.singleOccurrencesCount = res;
                return cqr;
            }
            
            
            inline CountQueriesResults countQueriesFor(const string& kmer) {
                uint_read_len_std kmerlen = kmer.length();
                if ((flatCqCache == 0) || (kmerlen > cqCacheMaxKLength))
                    return unknownCacheCount;

                uint_max idx = 0;
                for (uint_read_len_std k = 0; k < kmerlen; k++) {
                    int cqSymbol = symbolOrder[(unsigned char) kmer[k]];
                    if (cqSymbol == -1)
                        return unknownCacheCount;
                    idx = idx * cqSymbolsCount + cqSymbol;
                }

                if (kmerlen >= byteCacheKLength)
                    return convert(cqByteCache[kmerlen][idx]);
                else if (kmerlen >= shortCacheKLength)
                    return convert(cqShortCache[kmerlen][idx]);
                else
                    return cqCache[kmerlen][idx];
            }
            
        public:

            DefaultCountQueriesCache(ReadsSetProperties* readsSetProperties, uint_read_len_std cqCacheMaxKLength, uint_read_len_std shortCacheKLength, uint_read_len_std byteCacheKLength)
            {
                int j = 0;
                for(int i = 0; i < readsSetProperties->symbolsCount; i++)
                    if (readsSetProperties->symbolsList[i] != 'N')
                        symbolsList[j++] = readsSetProperties->symbolsList[i];
                symbolsList[j] = 0;
                
                cqSymbolsCount = j;
                generateSymbolOrder();
       
                this->cqCacheMaxKLength = cqCacheMaxKLength;
                this->shortCacheKLength = shortCacheKLength;
                this->byteCacheKLength = byteCacheKLength;
                
                flatCqCache = new CountQueriesResults[sizeOfFlatCqCache()]();
                flatCqShortCache = new unsigned short[sizeOfFlatCqShortCache()]();
                flatCqByteCache = new unsigned char[sizeOfFlatCqByteCache()]();
                
                buildArrays();
            }
            
            DefaultCountQueriesCache(std::istream& src) {
                int srchelper;
                src >> srchelper;
                cqSymbolsCount = srchelper;
                src.get(); // '/n'
                
                src >> symbolsList;
                src.get(); // '/n'
                
                src >> srchelper;
                cqCacheMaxKLength = srchelper;
                src.get(); // '/n'
                src >> srchelper;
                shortCacheKLength = srchelper;
                src.get(); // '/n'
                src >> srchelper;
                byteCacheKLength = srchelper;
                src.get(); // '/n'  
                
                flatCqCache = (CountQueriesResults*) PgSAHelpers::readArray(src, sizeOfFlatCqCache() * sizeof(CountQueriesResults));
                src.get(); // '/n'  
                flatCqShortCache = (unsigned short*) PgSAHelpers::readArray(src, sizeOfFlatCqShortCache() * sizeof(unsigned short));
                src.get(); // '/n'  
                flatCqByteCache = (unsigned char*) PgSAHelpers::readArray(src, sizeOfFlatCqByteCache() * sizeof(unsigned char));
                src.get(); // '/n'  
                
                generateSymbolOrder();
                buildArrays();
            }
            
            ~DefaultCountQueriesCache() {
                delete[]cqCache;
                delete[]cqShortCache;
                delete[]cqByteCache;
                delete[]flatCqCache;
                delete[]flatCqShortCache;
                delete[]flatCqByteCache;
            }
            
            string getTypeID() { return CACHETYPE_DEFAULT; };

            void write(std::ostream& dest) {
                dest << (int) cqSymbolsCount << "\n";
                dest << symbolsList << "\n";
                generateSymbolOrder();
                
                dest << (int) cqCacheMaxKLength << "\n";
                dest << (int) shortCacheKLength << "\n";
                dest << (int) byteCacheKLength << "\n";
                
                PgSAHelpers::writeArray(dest, flatCqCache, sizeOfFlatCqCache() * sizeof(CountQueriesResults));
                dest << "\n";
                PgSAHelpers::writeArray(dest, flatCqShortCache, sizeOfFlatCqShortCache() * sizeof(unsigned short));
                dest << "\n";
                PgSAHelpers::writeArray(dest, flatCqByteCache, sizeOfFlatCqByteCache() * sizeof(unsigned char));
                dest << "\n";
            }
            
            // Q2 - In how many reads does kmer occur?
            inline api_uint_reads_cnt countReadsImpl(const string& kmer) { 
                 return countQueriesFor(kmer).readsCount;
            };

            // Q4 - What is the number of occurrences of kmer?
            inline uint_max countOccurrencesImpl(const string& kmer) { 
                return countQueriesFor(kmer).occurrencesCount;
            };

            // Q6 - In how many reads does kmer occur only once?
            inline api_uint_reads_cnt countSingleOccurrencesImpl(const string& kmer) { 
                return countQueriesFor(kmer).singleOccurrencesCount;
            };
            
            inline const string getDescriptionImpl() { 
                string desc;
                uint_max size = sizeof(CountQueriesResults) * sizeOfFlatCqCache() +
                                sizeof(unsigned short) * sizeOfFlatCqShortCache() +
                                sizeof(unsigned char) * sizeOfFlatCqByteCache();
                desc = "Count Q CACHE " + toMB(size, 2) + " MB\tmax k: " +
                        toString(cqCacheMaxKLength) + 
                        " (16bit: " + toString(shortCacheKLength) + "-" + toString(byteCacheKLength - 1) + 
                        " / 8bit: " + toString(byteCacheKLength) + "-" + toString(cqCacheMaxKLength) + ")\n";                        
                return desc;
            };
    };
    
    typedef DefaultCountQueriesCache<uint_reads_cnt_std> DefaultCountQueriesCacheStandard;
    
    class GeneratedCountQueriesCache: public DefaultCountQueriesCacheStandard {
        
        private:
            
            void updateByteCache(uint_read_len_std kmerlen, uint_max idx, uint_reads_cnt_std addInHowManyReadsDoesFOccur,
                uint_max moreWhatIsTheNumberOccurrencesOfF, uint_reads_cnt_std lessInHowManyReadsDoesFOccurOnce) {
                if ((cqByteCache[kmerlen][idx] + addInHowManyReadsDoesFOccur >= exceedsUCharCache)
                    || (moreWhatIsTheNumberOccurrencesOfF + lessInHowManyReadsDoesFOccurOnce > 0)) {
                    cqByteCache[kmerlen][idx] = exceedsUCharCache;
                    return;
                }
                cqByteCache[kmerlen][idx] += addInHowManyReadsDoesFOccur;
            }

            inline void updateShortCache(uint_read_len_std kmerlen, uint_max idx, uint_reads_cnt_std addInHowManyReadsDoesFOccur,
                uint_max moreWhatIsTheNumberOccurrencesOfF, uint_reads_cnt_std lessInHowManyReadsDoesFOccurOnce) {
                if ((cqShortCache[kmerlen][idx] + addInHowManyReadsDoesFOccur >= exceedsUShortCache)
                    || (moreWhatIsTheNumberOccurrencesOfF + lessInHowManyReadsDoesFOccurOnce > 0)) {
                    cqShortCache[kmerlen][idx] = exceedsUShortCache;
                    return;
                }
                cqShortCache[kmerlen][idx] += addInHowManyReadsDoesFOccur;
            }

            inline void updateCache(uint_read_len_std kmerlen, uint_max idx, uint_reads_cnt_std addInHowManyReadsDoesFOccur,
                uint_max moreWhatIsTheNumberOccurrencesOfF, uint_reads_cnt_std lessInHowManyReadsDoesFOccurOnce) {
                if (kmerlen >= byteCacheKLength)
                    updateByteCache(kmerlen, idx, addInHowManyReadsDoesFOccur, moreWhatIsTheNumberOccurrencesOfF, lessInHowManyReadsDoesFOccurOnce);
                else if (kmerlen >= shortCacheKLength)
                    updateShortCache(kmerlen, idx, addInHowManyReadsDoesFOccur, moreWhatIsTheNumberOccurrencesOfF, lessInHowManyReadsDoesFOccurOnce);
                else {
                    cqCache[kmerlen][idx].readsCount += addInHowManyReadsDoesFOccur;
                    cqCache[kmerlen][idx].occurrencesCount += addInHowManyReadsDoesFOccur + moreWhatIsTheNumberOccurrencesOfF;
                    cqCache[kmerlen][idx].singleOccurrencesCount += addInHowManyReadsDoesFOccur - lessInHowManyReadsDoesFOccurOnce;
                }
            }
        
        public:
            
            template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
            GeneratedCountQueriesCache(DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>* pg, 
                        ReadsSetProperties* readsSetProperties,
                        uint_read_len_std cqCacheMaxKLength, uint_read_len_std shortCacheKLength, uint_read_len_std byteCacheKLength)
            : DefaultCountQueriesCache(readsSetProperties, cqCacheMaxKLength, shortCacheKLength, byteCacheKLength)
            { 
                clock_checkpoint();
                
                typedef typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type ReadsListClass;
                typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> ReadsList;
                ReadsList& readsList = *(pg->getReadsList());
                
                uint_max combinedCqCacheSize = sizeOfFlatCqCache(1, cqCacheMaxKLength);

                int_max* lastStart = new int_max[combinedCqCacheSize];
                int_max* lastEnd = new int_max[combinedCqCacheSize];

                for(uint_max i = 0; i < combinedCqCacheSize; i++) {
                    lastStart[i] = -1;
                    lastEnd[i] = -1;
                }

                uint_max cqCacheIdx[cqCacheMaxKLength + 1];
                cqCacheIdx[1] = 0;
                for (uint_read_len i = 2; i <= cqCacheMaxKLength; i++)
                    cqCacheIdx[i] = cqCacheIdx[i-1] + sizeOfCqCacheFor(i-1);

                int_max j_start = 0;
                int_max j_end = 0;
                for(uint_pg_len i = 0; i < pg->getLength(); i++) {
                    while (readsList.getReadPosition(j_start) + readsList.getMaxReadLength() <= i)
                        j_start++;
                    while (readsList.getReadPosition(j_end + 1) <= i)
                        j_end++;
                    int_max j = j_start;
                    uint_max idx = 0;
                    const char* pgenPtr = pg->getSuffix(i);
                    for (uint_read_len k = 1; k <= cqCacheMaxKLength; k++) {
                        if (readsList.getReadPosition(j_end) + readsList.getMaxReadLength() < i + k)
                            break;
                        int cqSymbol = symbolOrder[(unsigned char) *pgenPtr++];
                        if (cqSymbol == -1)
                            break;
                        idx = idx * cqSymbolsCount + cqSymbol;

                        uint_max flatIdx = cqCacheIdx[k] + idx;
                        if (lastEnd[flatIdx] >= j) {
                            if (lastStart[flatIdx] > j)
                                updateCache(k, idx, (j_end - lastEnd[flatIdx]), (lastEnd[flatIdx] - j) + 1,
                                             (lastEnd[flatIdx] - lastStart[flatIdx]) + 1);
                            else
                                updateCache(k, idx, (j_end - lastEnd[flatIdx]), (lastEnd[flatIdx] - j) + 1,
                                             (lastEnd[flatIdx] - j) + 1);
                            lastStart[flatIdx] = lastEnd[flatIdx] + 1;
                        } else {
                            updateCache(k, idx, (j_end - j) + 1, 0, 0);
                            lastStart[flatIdx] = j;
                        }

                        lastEnd[flatIdx] = j_end;

                        while (readsList.getReadPosition(j) + readsList.getMaxReadLength() <= i + k)
                            j++;
                    }
                }
                delete[]lastStart;
                delete[]lastEnd;
                
                cout << "Count Queries Cache done in " << clock_millis() << " msec\n\n";
            }
            
            static const uint_max bigPgLength = 300000000;
            
            template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
            GeneratedCountQueriesCache(DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>* pg, 
                        ReadsSetProperties* readsSetProperties)
            : GeneratedCountQueriesCache(pg, readsSetProperties, 13, pg->getLength()>bigPgLength?12:11, 13)
            {
            }
        
    };
}





#endif	/* DEFAULTCOUNTQUERIESCACHE_H */

