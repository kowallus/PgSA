#include "GreedySwipingDefaultOverlapPseudoGenomeGenerator.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::GreedySwipingDefaultOverlapPseudoGenomeGenerator(DefaultReadsSet* orgReadsSet)
    {
        if (!orgReadsSet->isReadLengthConstant())
            cout << "Unsupported: variable length reads :(";

        // warning: now generator handles orgReadsSet destruction
        this->defaultReadsSet = orgReadsSet;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::~GreedySwipingDefaultOverlapPseudoGenomeGenerator() {
        if (this->defaultReadsSet)
            delete(this->defaultReadsSet);
    }
    
        template<typename uint_read_len, typename uint_reads_cnt>
    string GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::getReadUpToOverlap(uint_reads_cnt incIdx) {
        return defaultReadsSet->getRead(incIdx - 1);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    uint_read_len GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::readLength(uint_reads_cnt incIdx) {
        return defaultReadsSet->readLength(incIdx - 1);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::readsTotal() {
        return defaultReadsSet->readsCount();
    }  
    template<typename uint_read_len, typename uint_reads_cnt>
    ReadsSetProperties* GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::getReadsSetProperties() {
        return defaultReadsSet->getReadsSetProperties();
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::compareReads(uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx) {
        return strcmp(defaultReadsSet->getRead(lIncIdx - 1).c_str(), defaultReadsSet->getRead(rIncIdx - 1).c_str());
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uchar GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::getSymbolOrderFromRead(uint_reads_cnt incIdx, uint_read_len offset) {
        return (uchar) getReadsSetProperties()->symbolOrder[(uchar) (defaultReadsSet->getRead(incIdx - 1).data()[offset])];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::compareSuffixes(uchar lSymOrder, uchar rSymOrder, uint_read_len offset) {
        uint_reads_cnt lIncIdx = sortedSuffixIdxs[ssiSymbolIdx[lSymOrder]];
        uint_reads_cnt rIncIdx = sortedSuffixIdxs[ssiSymbolIdx[rSymOrder]];
        return strcmp(defaultReadsSet->getRead(lIncIdx - 1).c_str() + offset, defaultReadsSet->getRead(rIncIdx - 1).c_str() + offset);
            
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::compareSuffixWithPrefix(uint_reads_cnt sufIncIdx, uint_reads_cnt preIncIdx, uint_read_len sufOffset) {
        return readsSufPreCmp(defaultReadsSet->getRead(sufIncIdx - 1).c_str() + sufOffset, defaultReadsSet->getRead(preIncIdx - 1).c_str());
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::updateSuffixQueue(uchar suffixGroup, uint_read_len suffixOffset) {
        if (ssiSymbolIdx[suffixGroup] < ssiSymbolEnd[suffixGroup]) {
            list<uchar>::reverse_iterator it = ssiOrder.rbegin();
            while (true) {
                if (it == ssiOrder.rend() || compareSuffixes(suffixGroup, *it, suffixOffset) >= 0) {
                    ssiOrder.insert(it.base(), suffixGroup);
                    break;
                }
                it++;
            }
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::initAndFindDuplicates() {

        for(uint_reads_cnt i = 1; i <= defaultReadsSet->readsCount(); i++)
            sortedReadsIdxs.push_back(i);
        
        PackedReadsComparator comparePacked = PackedReadsComparator(this);
        std::sort(sortedReadsIdxs.begin(), sortedReadsIdxs.end(), comparePacked);
        
        vector<uint_reads_cnt> sortedReadsLeft;
                
        sortedReadsLeft.push_back(*(sortedReadsIdxs.begin()));
        uchar curSymOrder = 0;
        
        for(typename vector<uint_reads_cnt>::iterator srIt = sortedReadsIdxs.begin(); srIt != sortedReadsIdxs.end();) {
            srIt++;
            if (srIt != sortedReadsIdxs.end() && compareReads(*(srIt-1),*srIt) == 0) 
                this->setReadSuccessor(*(srIt-1), *srIt, defaultReadsSet->maxReadLength());
            else {
                sortedSuffixIdxs.push_back(*(srIt-1));
                uchar firstSymbolOrder = getSymbolOrderFromRead(*(srIt-1), 0);
                if (curSymOrder != firstSymbolOrder) {
                    uint_reads_cnt ssiIndex = sortedSuffixIdxs.size() - 1;
                    ssiSymbolEnd[curSymOrder] = ssiIndex;
                    ssiSymbolIdx[firstSymbolOrder] = ssiIndex;
                    curSymOrder = firstSymbolOrder;
                }
                if (srIt != sortedReadsIdxs.end())
                    sortedReadsLeft.push_back(*srIt);
            }
        }  
        ssiSymbolEnd[curSymOrder] = sortedSuffixIdxs.size();
        sortedReadsIdxs = sortedReadsLeft;
        
        cout << "Found " << (readsTotal() - this->readsLeft) << " duplicates.\n";
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>::findOverlappingReads() {

        initAndFindDuplicates();
        cout << "Start overlapping.\n";
    
        for(int i = 1 ; i < defaultReadsSet->maxReadLength(); i++) {
            
            vector<uint_reads_cnt> sortedReadsLeft;
            sortedReadsLeft.reserve(this->readsLeft);
            vector<uint_reads_cnt> sortedSuffixLeft;
            sortedSuffixLeft.reserve(this->readsLeft);
            vector<uint_reads_cnt> ssiSymbolIdxLeft(UCHAR_MAX, 0);
            vector<uint_reads_cnt> ssiSymbolEndLeft(UCHAR_MAX, 0);
            uchar curSymOrder = 0;
            
            for(int j = 0; j < getReadsSetProperties()->symbolsCount; j++) 
                updateSuffixQueue(j, i);
                  
            typename vector<uint_reads_cnt>::iterator preIt = sortedReadsIdxs.begin();
            while (!ssiOrder.empty() || (preIt != sortedReadsIdxs.end())) {
                if (ssiOrder.empty()) 
                    sortedReadsLeft.push_back(*(preIt++));
                else {
                    uchar j = ssiOrder.front();
                    uint_reads_cnt sufIdx = sortedSuffixIdxs[ssiSymbolIdx[j]];
                              
                    if (preIt != sortedReadsIdxs.end()) {
                        int cmpRes = -1;
                        
                        typename vector<uint_reads_cnt>::iterator curPreIt = preIt;
                        while (preIt != sortedReadsIdxs.end()) {
                            if ((cmpRes = compareSuffixWithPrefix(sufIdx, *preIt, i)) != 0)
                                break;
                            if ((sufIdx != *preIt) && (!this->isHeadOf(sufIdx, *preIt)))
                                break;
                            cmpRes = -1;
                            preIt++;
                        }
                    
                        if (cmpRes)
                            preIt = curPreIt;
                        else {
                            uint_reads_cnt preIdx = *preIt;
                            while (preIt > curPreIt) {
                                *preIt = *(preIt - 1);
                                preIt--;
                            }
                            *preIt = preIdx;
                        }
                            
                        if (cmpRes == 0) 
                            this->setReadSuccessor(sufIdx, *(preIt++), defaultReadsSet->maxReadLength() - i);
                        else if (cmpRes > 0) {
                            sortedReadsLeft.push_back(*(preIt++));
                            continue;
                        } else {
                            sortedSuffixLeft.push_back(sufIdx);
                            uchar ithSymbolOrder = getSymbolOrderFromRead(sufIdx, i);
                            if (curSymOrder != ithSymbolOrder) {
                                uint_reads_cnt ssiIndex = sortedSuffixLeft.size() - 1;
                                ssiSymbolEndLeft[curSymOrder] = ssiIndex;
                                ssiSymbolIdxLeft[ithSymbolOrder] = ssiIndex;
                                curSymOrder = ithSymbolOrder;
                            }
                        }
                    }
                    
                    ssiOrder.pop_front();
                    ssiSymbolIdx[j]++;
                    updateSuffixQueue(j, i);
                }
                               
            }
            
            ssiSymbolEndLeft[curSymOrder] = sortedSuffixLeft.size();
            sortedReadsIdxs = sortedReadsLeft;
            sortedSuffixIdxs = sortedSuffixLeft;
            ssiSymbolIdx = ssiSymbolIdxLeft;
            ssiSymbolEnd = ssiSymbolEndLeft;
        
            cout << this->readsLeft << " reads left after " << (uint_read_len_max) (defaultReadsSet->maxReadLength() - i) << " overlap\n";
        }
        
        cout << this->countComponents() << " pseudo-genome components\n";
        
        cout << "Overlapping done in " << clock_millis() << " msec\n\n";     
    }
    
// FACTORY

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeGeneratorBase* GreedySwipingDefaultOverlapPseudoGenomeGeneratorFactory::getGeneratorFullTemplate(DefaultReadsSet* readsSet) {
        return new GreedySwipingDefaultOverlapPseudoGenomeGenerator<uint_read_len, uint_reads_cnt>(readsSet);
    }

    template<typename uint_read_len>
    PseudoGenomeGeneratorBase* GreedySwipingDefaultOverlapPseudoGenomeGeneratorFactory::getGeneratorPartialTemplate(DefaultReadsSet* readsSet) {

        if (isReadsCountStd(readsSet->readsCount()))
            return getGeneratorFullTemplate<uint_read_len, uint_reads_cnt_std>(readsSet);
        else
            cout << "UNSUPPORTED READS COUNT!!!???";

        return 0;
    }

    PseudoGenomeGeneratorBase* GreedySwipingDefaultOverlapPseudoGenomeGeneratorFactory::getGenerator(string readsFile, string pairFile) {

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
