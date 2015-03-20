#include "GreedyVerticalOverlapPseudoGenomeGenerator.h"
#include "../../readsset/DefaultReadsSet.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::GreedyVerticalOverlapGeneratorTemplate(DefaultReadsSet* orgReadsSet)
    {
        if (!orgReadsSet->isReadLengthConstant())
            cout << "Unsupported: variable length reads :(";

        // warning: now generator handles orgReadsSet destruction
        this->orgReadsSet = orgReadsSet;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::~GreedyVerticalOverlapGeneratorTemplate() {
         delete(this->orgReadsSet);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    string GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getReadUpToOverlap(uint_reads_cnt incIdx) {
        return orgReadsSet->getRead(incIdx - 1);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    uint_read_len GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::readLength(uint_reads_cnt incIdx) {
        return orgReadsSet->readLength(incIdx - 1);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::readsTotal() {
        return orgReadsSet->readsCount();
    }  
    template<typename uint_read_len, typename uint_reads_cnt>
    ReadsSetProperties* GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getReadsSetProperties() {
        return orgReadsSet->getReadsSetProperties();
    }

    
    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::findOverlappingReads() {

        populateReadsSet();
        
        cout << "Start overlapping.\n";
        
        overlapMultisetVertical(this->orgReadsSet->maxReadLength(), 1);

        readsSet.clear();

        cout << "Overlapping done in " << clock_millis() << " msec\n\n";     
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::populateReadsSet() {
        typename std::multiset<ReadWithIdx>::iterator it;
        uint_reads_cnt i = 0;
        while (readsSet.empty() && ++i <= this->readsTotal())
            if (!this->hasPredecessor(i))
                it = readsSet.insert({i, getReadUpToOverlap(i).c_str()});
        while (++i <= readsTotal())
            if (!this->hasPredecessor(i))
                it = readsSet.insert(it, {i, getReadUpToOverlap(i).c_str()});

        cout << "Put " << readsSet.size() << " reads in multiset.\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::matchInReadsSet(uint_reads_cnt suffixIdx, uint_read_len overlap) {
        const char* suffix = getReadUpToOverlap(suffixIdx).c_str() + (readLength(suffixIdx) - overlap);

        typename std::multiset<ReadWithIdx>::iterator it = (readsSet.lower_bound({0, suffix}));
        uint_reads_cnt prefixReadIdx;
        bool match = false;
        while(it != readsSet.end()
                && (prefixReadIdx = (*it).incIdx)
                && (match = (readsSufPreCmp(suffix, (*it).read) == 0))
                && ((suffixIdx == prefixReadIdx) || this->isHeadOf(suffixIdx, prefixReadIdx)))
            it++;
        if (match && (it != readsSet.end())) {
            this->setReadSuccessor(suffixIdx, prefixReadIdx, overlap);

            readsSet.erase(it);
            return prefixReadIdx;
        }
        return 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedyVerticalOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::overlapMultisetVertical(uint_read_len maxOverlap, uint_read_len overlapLimit) {

        for (uint_read_len j = maxOverlap; j >= overlapLimit; j--) {
            if (!this->readsLeft)
                break;

            // find biggest overlap
            for (uint_reads_cnt i = 1; i <= readsTotal(); i++) {
                if (this->hasSuccessor(i))
                    continue;
                matchInReadsSet(i, j);
            }
            cout << this->readsLeft << " reads left after " << (uint_read_len_max) j << " overlap\n";
        }

//        cout << this->readsLeft << " reads left after " << (uint_read_len_max) overlapLimit << " overlap\n";
        cout << this->countComponents() << " pseudo-genome components\n";
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
