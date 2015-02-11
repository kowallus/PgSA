#ifndef LISTOFCONSTANTLENGTHREADS_H_INCLUDED
#define LISTOFCONSTANTLENGTHREADS_H_INCLUDED

#define ORIGINAL_INDEX_OFFSET (FLAGS_OFFSET+1)

#include "ReadsListInterface.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {
    
    template <  typename uint_read_len,
                typename uint_reads_cnt,
                typename uint_pg_len,
                unsigned char LIST_ELEMENT_SIZE,              // size of list element in uchar
                uchar FLAGS_OFFSET          // offset of flags in uchar (flags+1 is original read index offset)
                >
    class ListOfConstantLengthReads: public GeneratedReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len,
                                                ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>>
                                    ,public ReadsListIteratorInterface<uint_read_len, uint_reads_cnt, 
                                                ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>>
    {
        private:
            
            typedef ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET> ThisReadsListType;
            typedef ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ThisReadsListType> ReadsList;
            
            static int elementsCompare (const void * a, const void * b) {
                if (*(uint_pg_len*)a > *(uint_pg_len*)b)
                    return 1;
                if (*(uint_pg_len*)a < *(uint_pg_len*)b)
                    return -1;
                return 0;
            }

            uchar* pgReadsList = 0;
            uchar* pgReadsListEnd = 0;

            const uint_read_len readLength;

            uint_reads_cnt readsCount;
            uint_pg_len pseudoGenomeLength;

            uint_read_len duplicateFilterK;
            
            uint_reads_cnt* readsListIdx = 0; // conversion from original index to reads list index
            
            ///////////////////////////
            // GENERATION VARIABLES
            
            uint_reads_cnt curRawIdx;
            uint_pg_len maxPos;
            bool isSortRequired;

            void generateReverseIndex() {
                readsListIdx = new uint_reads_cnt[readsCount];
                for (uint_reads_cnt i = 0; i < readsCount; i++)
                    readsListIdx[this->getReadOriginalIndexImpl(i)] = i;                
            }
            
            ///////////////////////////
            // ITERATION VARIABLES
                                    
            uint_pg_len suffixPos;
            int_max guard; // negative possible
            
            uint_read_len guardOffset;
            bool skipDuplicateFilter;
            
            uchar* itAddress;
            
            vector<uchar*> occurFlagsReadsList;
            vector<pair<uchar*, uint_read_len>> occurFlagsOccurrencesList;
            
            ////////////////////////////
            // HELPER ROUTINES
            
            inline uchar* idx2Address(uint_reads_cnt idx) {
                return pgReadsList + ((uint_max) idx) * LIST_ELEMENT_SIZE ;
            };
            
            inline uchar& flagsVariableByAddress(uchar* address) { return ( *(address + FLAGS_OFFSET)); };
            
            inline uint_pg_len getReadPositionByAddress(uchar* address) {
                // casting address to larger int (uint_max*) improved caching
                return (uint_pg_len) *((uint_max*) address);
                //return *((uint_pg_len*) address);
            };
            
            inline uint_reads_cnt getReadOriginalIndexByAddress(uchar* address) {
                return (*((uint_reads_cnt*) (address + ORIGINAL_INDEX_OFFSET)));
            };
            
            inline bool hasDuplicateFilterFlagByAddress(uchar* address) {
                return flagsVariableByAddress(address) & DUPINREADS_FLAG;
            };
            
            inline void setDuplicateFilterFlagByAddress(uchar* address) { 
                flagsVariableByAddress(address) |= DUPINREADS_FLAG;
            };

            inline bool hasOccurFlagByAddress(uchar* address) {
                return flagsVariableByAddress(address) & OCCUR_FLAG;
            }

            inline bool hasOccurOnceFlagByAddress(uchar* address) {
                return (flagsVariableByAddress(address) & OCCURFLAGS_MASK) == OCCUR_FLAG;
            }

            inline void setOccurFlagByAddress(uchar* address) {
                flagsVariableByAddress(address) |= OCCUR_FLAG;
            }

            inline void clearOccurFlagsByAddress(uchar* address) {
                flagsVariableByAddress(address) &= ~(OCCURFLAGS_MASK);
            }

            inline void setOccurOnceFlagByAddress(uchar* address) {
                flagsVariableByAddress(address) |= (((flagsVariableByAddress(address) & OCCUR_FLAG) << 1) | OCCUR_FLAG);
            }
            
            inline uchar& flagsVariable(uint_reads_cnt idx) { return ( *(idx2Address(idx) + FLAGS_OFFSET)); };
            
        public:

            ListOfConstantLengthReads(uint_read_len readLength, uint_reads_cnt readsCount, uint_pg_len pseudoGenomeLength)
            : readLength(readLength)
            {
                this->readsCount = readsCount;
                this->pseudoGenomeLength = pseudoGenomeLength;
                this->pgReadsList = new uchar[( readsCount + 1) * (uint_max) LIST_ELEMENT_SIZE];
                this->curRawIdx = 0;
                this->maxPos = 0;
                this->isSortRequired = false;
                this->duplicateFilterK = -1;
            }

            ListOfConstantLengthReads(uint_read_len readLength, std::istream& src)
            : readLength(readLength) {
                src >> readsCount;
                src >> pseudoGenomeLength;
                int srchelper;
                src >> srchelper;
                duplicateFilterK = srchelper;
                src.get(); //"\n";
                this->pgReadsList = (uchar*) PgSAHelpers::readArray(src, sizeof(uchar) * (readsCount + 1) * LIST_ELEMENT_SIZE);

                this->isSortRequired = false;
                generateReverseIndex();
                this->pgReadsListEnd = this->pgReadsList + (uint_max) readsCount * LIST_ELEMENT_SIZE;
            }

            void writeImpl(std::ostream& dest) {
                if (isSortRequired)
                    cout << "WARNING: Reads list invalid!";
                dest << readsCount << "\n";
                dest << pseudoGenomeLength << "\n";
                dest << (int) duplicateFilterK << "\n";
                PgSAHelpers::writeArray(dest, this->pgReadsList, sizeof(uchar) * (readsCount + 1) * LIST_ELEMENT_SIZE );
            };

            virtual ~ListOfConstantLengthReads() {
                delete[] (this->pgReadsList);
                if (this->readsListIdx != 0)
                    delete[] (this->readsListIdx);
            };
            
            inline const uchar getListElementSizeImpl() { return LIST_ELEMENT_SIZE; };

            inline uint_pg_len getReadPositionImpl(uint_reads_cnt idx) { return getReadPositionByAddress(idx2Address(idx)); };

            inline uint_reads_cnt getReadsCountImpl() { return readsCount; };
            
            inline uint_read_len getMaxReadLengthImpl() { return readLength; };
            inline uint_read_len getReadLengthImpl(uint_reads_cnt) { return readLength; };

            inline uint_reads_cnt getReadOriginalIndexImpl(uint_reads_cnt idx) { return getReadOriginalIndexByAddress(idx2Address(idx)); };

            inline uint_reads_cnt getReadsListIndexOfOriginalIndexImpl(uint_reads_cnt originalIdx) { 
                return readsListIdx[originalIdx];
            }
            
            inline uint_read_len getDuplicateFilterKmerLengthImpl() { return duplicateFilterK; };
            inline void setDuplicateFilterKmerLengthImpl(uint_read_len kLength) { duplicateFilterK = kLength; };

            inline bool hasDuplicateFilterFlagImpl(uint_reads_cnt idx) { return hasDuplicateFilterFlagByAddress(idx2Address(idx)); };

            inline void setDuplicateFilterFlagImpl(uint_reads_cnt idx) { setDuplicateFilterFlagByAddress(idx2Address(idx)); };

            /////////////////////////////////////////
            // OCCUR FLAGS MANAGEMENT
            
            inline bool hasOccurFlagImpl(uint_reads_cnt idx) { return hasOccurFlagByAddress(idx2Address(idx)); };

            inline bool hasOccurOnceFlagImpl(uint_reads_cnt idx) { return hasOccurOnceFlagByAddress(idx2Address(idx)); };

            inline void setOccurFlagImpl(uint_reads_cnt idx) { setOccurFlagByAddress(idx2Address(idx)); };

            inline void clearOccurFlagsImpl(uint_reads_cnt idx) { clearOccurFlagsByAddress(idx2Address(idx)); };

            inline void setOccurOnceFlagImpl(uint_reads_cnt idx) { setOccurOnceFlagByAddress(idx2Address(idx)); };

            ///////////////////////////
            // ITERATION ROUTINES
                                    
            inline void initIterationImpl(const uint_read_len& kmerLength) { 
                guardOffset = readLength - kmerLength;
                skipDuplicateFilter = kmerLength < duplicateFilterK;
                itAddress = pgReadsList;
            };
            
            inline void setIterationPositionImpl(const RPGOffset<uint_read_len, uint_reads_cnt>& saPos) {
                itAddress = idx2Address(saPos.readListIndex);
                suffixPos = saPos.offset + this->getReadPositionByAddress(itAddress); 
                guard = (int_max) suffixPos - guardOffset;
                itAddress += LIST_ELEMENT_SIZE;;
            };

            inline void setIterationPositionImpl(const RPGOffset<uint_read_len, uint_reads_cnt>& saPos, uchar shift) {
                itAddress = idx2Address(saPos.readListIndex);
                suffixPos = saPos.offset + this->getReadPositionByAddress(itAddress) - shift;
                guard = (int_max) suffixPos - guardOffset;

                while ((itAddress >= pgReadsList) && (suffixPos < getReadPositionByAddress(itAddress)))
                    itAddress -= LIST_ELEMENT_SIZE;
                itAddress += LIST_ELEMENT_SIZE;
            };
            
            inline bool moveNextImpl() {
                return ((itAddress -= LIST_ELEMENT_SIZE) >= pgReadsList) && ((int_max) this->getReadPositionByAddress(itAddress) >= guard);
            };
            
             inline uint_reads_cnt getReadOriginalIndexImpl() {
                return getReadOriginalIndexByAddress(itAddress);
            };
            
            inline uint_read_len getOccurrenceOffsetImpl() { 
                return suffixPos - getReadPositionByAddress(itAddress);
            };
            
            inline uint_flatten_occurrence_max getFlattenOccurrenceImpl() { 
                return ((uint_flatten_occurrence_max) this->getReadOriginalIndexImpl()) * this->getMaxReadLength() + this->getOccurrenceOffsetImpl();
            };
            
            inline bool hasDuplicateFilterFlagImpl() { 
                return (skipDuplicateFilter || hasDuplicateFilterFlagByAddress(itAddress));
            };

            inline bool hasOccurFlagImpl() { 
                return hasOccurFlagByAddress(itAddress);
            };
            
            inline bool hasOccurOnceFlagImpl() { 
                return hasOccurOnceFlagByAddress(itAddress);
            };

            inline void setOccurFlagImpl() { 
                setOccurFlagByAddress(itAddress);
                occurFlagsReadsList.push_back(itAddress);
            };

            inline uint_reads_cnt clearAllOccurFlagsImpl() { 
                uint_reads_cnt readsWithOccurFlagCount = occurFlagsReadsList.size();
                for(typename vector<uchar*>::iterator it = occurFlagsReadsList.begin(); it != occurFlagsReadsList.end(); ++it) 
                    clearOccurFlagsByAddress(*it);
                
                occurFlagsReadsList.clear();
                
                return readsWithOccurFlagCount;
            };

            inline void setOccurOnceFlagAndPushReadImpl() { 
                setOccurOnceFlagByAddress(itAddress);
                occurFlagsReadsList.push_back(itAddress);
            };
            
            inline void setOccurOnceFlagAndPushOccurrenceImpl() { 
                setOccurOnceFlagByAddress(itAddress);
                occurFlagsOccurrencesList.push_back( { itAddress, getOccurrenceOffsetImpl()} );
            };
            
            template<typename api_uint_reads_cnt>
            inline void clearAllOccurFlagsAndPushReadsWithSingleOccurrenceImpl(vector<api_uint_reads_cnt>& reads) { 
                for(typename vector<uchar*>::iterator it = occurFlagsReadsList.begin(); it != occurFlagsReadsList.end(); ++it) {
                    if (hasOccurOnceFlagByAddress(*it))
                        reads.push_back(getReadOriginalIndexByAddress(*it));
                    clearOccurFlagsByAddress(*it);
                };
                
                occurFlagsReadsList.clear();
            };
            
            inline uint_reads_cnt clearAllOccurFlagsAndCountSingleOccurrencesImpl() { 
                uint_reads_cnt count = 0;
                for(typename vector<uchar*>::iterator it = occurFlagsReadsList.begin(); it != occurFlagsReadsList.end(); ++it) {
                    if (hasOccurOnceFlagByAddress(*it))
                        count++;
                    clearOccurFlagsByAddress(*it);
                }
                
                occurFlagsReadsList.clear();
                
                return count;
            };
            
            template<typename api_uint_reads_cnt, typename api_uint_read_len>
            inline void clearAllOccurFlagsAndPushSingleOccurrencesImpl(vector<pair<api_uint_reads_cnt, api_uint_read_len>>& occurrences) { 
                for(typename vector<pair<uchar*, uint_read_len>>::iterator it = occurFlagsOccurrencesList.begin(); it != occurFlagsOccurrencesList.end(); ++it) {
                    if (this->hasOccurOnceFlagByAddress((*it).first))
                        occurrences.push_back({ getReadOriginalIndexByAddress((*it).first), (*it).second });
                    this->clearOccurFlagsByAddress((*it).first);
                }
                
                occurFlagsOccurrencesList.clear();
            };

            inline void clearAllOccurFlagsAndPushSingleOccurrencesFlattenImpl(vector<uint_flatten_occurrence_max>& flattenOccurrences) {
                for(typename vector<pair<uchar*, uint_read_len>>::iterator it = occurFlagsOccurrencesList.begin(); it != occurFlagsOccurrencesList.end(); ++it) {
                    if (this->hasOccurOnceFlagByAddress((*it).first))
                        flattenOccurrences.push_back( ((uint_flatten_occurrence_max) this->getReadOriginalIndexByAddress((*it).first)) * this->getMaxReadLength() + (*it).second );
                    this->clearOccurFlagsByAddress((*it).first);
                }
                
                occurFlagsOccurrencesList.clear();                        
            };            
            
            /////////////////////////////////
            // GENERATION ROUTINES
           
            inline void addImpl(uint_pg_len pos, uint_read_len len, uint_reads_cnt idx) {
                if (maxPos <= pos)
                    maxPos = pos;
                else
                    isSortRequired = true;
                *((uint_pg_len*) (this->pgReadsList + curRawIdx)) = pos;
                *((uint_reads_cnt*) (this->pgReadsList + curRawIdx + ORIGINAL_INDEX_OFFSET)) = idx;
                curRawIdx += LIST_ELEMENT_SIZE;
            }

            void validateImpl() {
                if (curRawIdx != readsCount * LIST_ELEMENT_SIZE)
                    cout << "WARNING: Reads list generation failed: " << curRawIdx / LIST_ELEMENT_SIZE << " elements instead of " << readsCount << "\n";

                if (isSortRequired) {
                    qsort(this->pgReadsList, readsCount, sizeof(uchar) * LIST_ELEMENT_SIZE, elementsCompare);
                    isSortRequired = false;
                }

                // add guard
                 *((uint_pg_len*) (this->pgReadsList + curRawIdx)) = this->pseudoGenomeLength;
            }

    };
    
    template <  typename uint_read_len,
             typename uint_reads_cnt,
             typename uint_pg_len,
             unsigned char LIST_ELEMENT_SIZE,              // size of list element in uchar
             uchar FLAGS_OFFSET          // offset of flags in uchar (flags+1 is original read index offset)
             > 
    struct ReadsListIteratorFactoryTemplate<ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>> {
        
        typedef ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET> ReadsListClass;
        typedef ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET> ReadsListIteratorClass;
        
        //FIXME: How to handle destruction? ... and get rid of *
        static ReadsListIteratorClass* getReadsListIterator(ReadsListClass& readsList) {
            return &readsList;
        };
    };
    
}
#endif // LISTOFCONSTANTLENGTHREADS_H_INCLUDED
