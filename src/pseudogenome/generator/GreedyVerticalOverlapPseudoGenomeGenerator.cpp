#include "GreedyVerticalOverlapPseudoGenomeGenerator.h"
#include "../../readsset/DefaultReadsSet.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::GreedyVerticalOverlapGeneratorTemplate(DefaultReadsSet* orgReadsSet) {
        if (!orgReadsSet->isReadLengthConstant())
            cout << "Unsupported: variable length reads :(";

        // warning: now generator handles orgReadsSet destruction
        this->orgReadsSet = orgReadsSet;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::~GreedyVerticalOverlapGeneratorTemplate() {
         delete(orgReadsSet);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::init() {   
        prevRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        nextRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        overlap = (uint_read_len*) calloc(readsTotal() + 1, sizeof(uint_read_len));
        headRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        readsLeft = readsTotal();
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::dispose() {   
        free(prevRead);
        free(nextRead);
        free(overlap);
        free(headRead);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeBase* GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::generatePseudoGenomeBase() {

        clock_checkpoint();
      
        init();  
        populateReadsSet();
        
        cout << "Start overlapping.\n";
        
        overlapMultisetVertical(orgReadsSet->maxReadLength(), 1);

        readsSet.clear();

        cout << "Overlapping done in " << clock_millis() << " msec\n\n";
        
        pseudoGenomeLength = countPseudoGenomeLength();
        quick_stats();

        PseudoGenomeBase* pgb = 0;
        
        if (isPseudoGenomeLengthStandardVirtual())
            pgb = assemblePseudoGenomeTemplate<uint_pg_len_std>();
        else if (isPseudoGenomeLengthMaximalVirtual())
            pgb = assemblePseudoGenomeTemplate<uint_pg_len_max>();
        else
            cout << "Unsupported: pseudo genome length :(";
        
        dispose();
        
        return pgb;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::populateReadsSet() {
        typename std::multiset<ReadWithIdx>::iterator it;
        uint_reads_cnt i = 0;
        while (readsSet.empty() && ++i <= readsTotal())
            if (!hasPredecessor(i))
                it = readsSet.insert({i, getRead(i).c_str()});
        while (++i <= readsTotal())
            if (!hasPredecessor(i))
                it = readsSet.insert(it, {i, getRead(i).c_str()});

        cout << "Put " << readsSet.size() << " reads in multiset.\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::matchInReadsSet(uint_reads_cnt suffixIdx, uint_read_len overlap) {
        const char* suffix = getRead(suffixIdx).c_str() + (readLength(suffixIdx) - overlap);

        typename std::multiset<ReadWithIdx>::iterator it = (readsSet.lower_bound({0, suffix}));
        uint_reads_cnt prefixReadIdx;
        bool match = false;
        while(it != readsSet.end()
                && (prefixReadIdx = (*it).incIdx)
                && (match = (readsSufPreCmp(suffix, (*it).read) == 0))
                && ((suffixIdx == prefixReadIdx) || isHeadOf(suffixIdx, prefixReadIdx)))
            it++;
        if (match && (it != readsSet.end())) {
            setReadSuccessor(suffixIdx, prefixReadIdx, overlap);

            readsSet.erase(it);
            return prefixReadIdx;
        }
        return 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::overlapMultisetVertical(uint_read_len maxOverlap, uint_read_len overlapLimit) {

        for (uint_read_len j = maxOverlap; j >= overlapLimit; j--) {
            if (!readsLeft)
                break;

            // find biggest overlap
            for (uint_reads_cnt i = 1; i <= readsTotal(); i++) {
                if (hasSuccessor(i))
                    continue;
                matchInReadsSet(i, j);
            }
//            cout << readsLeft << " reads left after " << (uint_read_len_max) j << " overlap\n";
        }

//        cout << readsLeft << " reads left after " << (uint_read_len_max) overlapLimit << " overlap\n";
        cout << countComponents() << " pseudo-genome components\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::setReadSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght) {
        nextRead[curIdx] = nextIdx;
        overlap[curIdx] = overlapLenght;
        if (headRead[nextIdx] == 0)
            headRead[curIdx] = nextIdx;
        else
            headRead[curIdx] = headRead[nextIdx];
        prevRead[nextIdx] = curIdx;
        readsLeft--;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getHead(uint_reads_cnt idx) {
        if (headRead[idx] == 0)
            return idx;

        uint_reads_cnt headIdx = idx;
        while(headRead[headIdx] != 0)
            headIdx = headRead[headIdx];
        headRead[idx] = headIdx;
        return headIdx;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_pg_len_max GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::countPseudoGenomeLength() {
        uint_pg_len_max len = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (overlap[i] < readLength(i))
                len += (readLength(i) - overlap[i]);

        return len;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::countComponents() {
        uint_reads_cnt count = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (!hasPredecessor(i) && hasSuccessor(i))
                count++;
        return count;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::countSingles() {
        uint_reads_cnt count = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (!hasPredecessor(i) && !hasSuccessor(i))
                count++;
        return count;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::quick_stats() {

        cout << pseudoGenomeLength << " bytes after overlapping\n";
        cout << countComponents() << " pseudo-genome components\n";
        cout << countSingles() << " single reads\n";
    }

// HELPERS

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead) {
        while (*suffixPart) {
            if (*suffixPart > *prefixRead)
                return 1;
            if (*suffixPart++ < *prefixRead++)
                return -1;
        }
        return 0;
    }

// FACTORY

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeGeneratorBase* GreedyVerticalOverlapPseudoGenomeGeneratorFactory::getGeneratorFullTemplate(DefaultReadsSet* readsSet) {
        return new GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>(readsSet);
    }

    template<typename uint_read_len>
    PseudoGenomeGeneratorBase* GreedyVerticalOverlapPseudoGenomeGeneratorFactory::getGeneratorPartialTemplate(DefaultReadsSet* readsSet) {

        if (isReadsCountStd(readsSet->readsCount()))
            return getGeneratorFullTemplate<uint_read_len, uint_reads_cnt_std>(readsSet);
        else
            cout << "UNSUPPORTED READS COUNT!!!???";

        return 0;
    }

    PseudoGenomeGeneratorBase* GreedyVerticalOverlapPseudoGenomeGeneratorFactory::getGenerator(string readsFile, string pairFile) {

        cout << "Reading reads set\n";
        
        // readsSet will be freed during generator destruction.        
        DefaultReadsSet* readsSet = DefaultReadsSet::readReadsSet(readsFile, pairFile);
          
        readsSet->printout();

        PseudoGenomeGeneratorBase* generatorBase = 0;

        if (isReadLengthMin(readsSet->maxReadLength()))
            generatorBase = getGeneratorPartialTemplate<uint_read_len_min>(readsSet);
        else if (isReadLengthStd(readsSet->maxReadLength()))
            generatorBase = getGeneratorPartialTemplate<uint_read_len_std>(readsSet);
        else cout << "UNSUPPORTED READS LENGTH!!!";
        
        return generatorBase;
    }
}
