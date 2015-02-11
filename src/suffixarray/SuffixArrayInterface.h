#ifndef SUFFIXARRAYINTERFACETEMPLATE_H_INCLUDED
#define SUFFIXARRAYINTERFACETEMPLATE_H_INCLUDED

#include "SuffixArrayBase.h"
#include "iterator/OccurrencesIteratorInterface.h"
#include "../readsset/ReadsSetInterface.h"

namespace PgSAIndex {

    typedef void* sa_pos_addr;

    template <class SuffixArrayClass> 
    struct SuffixArrayTraits{
        typedef int ReadsListClass; // wrong type !!!
        typedef int ReadsListIteratorClass; // wrong type !!!
        typedef int OccurrencesIteratorClass; // wrong type !!!
    };
    
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, 
                class SuffixArrayClass >
    class SuffixArrayInterface
    {        
        public:
               
            virtual ~SuffixArrayInterface() {};

            typedef OccurrencesIteratorInterface<uint_read_len, uint_reads_cnt, typename SuffixArrayTraits<SuffixArrayClass>::ReadsListIteratorClass, typename SuffixArrayTraits<SuffixArrayClass>::OccurrencesIteratorClass> OccurrencesIterator;
            
            inline OccurrencesIterator& getKmerOccurrencesIterator(const string& kmer) { return static_cast<SuffixArrayClass*>(this)->getKmerOccurrencesIteratorImpl(kmer); };
            inline OccurrencesIterator& getKmerOccurrencesIterator(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength) { return static_cast<SuffixArrayClass*>(this)->getKmerOccurrencesIteratorImpl(originalIdx, pos, kmerLength); };
            
            inline void findKmerRange(const char* kmer, const uint_read_len& kmerLength, SARange<uint_pg_len>& range) { static_cast<SuffixArrayClass*>(this)->findKmerRangeImpl(kmer, kmerLength, range); };

            inline const RPGOffset<uint_read_len, uint_reads_cnt> getPosition(const uint_pg_len posIdx) { return static_cast<SuffixArrayClass*>(this)->getPositionImpl(posIdx); };
            
            inline const string getSuffix(const uint_pg_len posIdx, const uint_pg_len length) { return static_cast<SuffixArrayClass*>(this)->getSuffixImpl(posIdx, length); };

            inline const string getDescription() { return static_cast<SuffixArrayClass*>(this)->getDescriptionImpl(); };
            
            ReadsSetInterface<uint_read_len, uint_reads_cnt>* getReadsSet() { return static_cast<SuffixArrayClass*>(this)->getReadsSetImpl(); };
            
            static SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>* castBase(SuffixArrayBase* base) {
                // TODO: validate
                return dynamic_cast<SuffixArrayInterface<uint_read_len, uint_reads_cnt, uint_pg_len, SuffixArrayClass>*>(base);
            }
    };

}

#endif // SUFFIXARRAYINTERFACETEMPLATE_H_INCLUDED
