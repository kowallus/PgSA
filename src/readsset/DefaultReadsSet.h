#ifndef DEFAULTREADSSET_H_INCLUDED
#define DEFAULTREADSSET_H_INCLUDED

#include <algorithm>
#include "ReadsSetBase.h"
#include "iterator/ReadsSetIterator.h"
#include "ReadsSetInterface.h"

namespace PgSAReadsSet {

    class DefaultReadsSet: public ReadsSetBase, public ReadsSetInterface< uint_read_len_max, uint_reads_cnt_max >
    {
        private:

            vector<string> reads;

        public:

            template<class ReadsSourceIterator>
            DefaultReadsSet(ReadsSourceIterator*);

            virtual ~DefaultReadsSet() {};

            inline uint_read_len_max maxReadLength() { return properties->maxReadLength; };
            inline uint_reads_cnt_max readsCount() { return properties->readsCount; };

            inline bool isReadLengthConstant() { return properties->constantReadLength; };

            inline const string getRead(uint_reads_cnt_max i) { return reads[i];};
            inline uint_read_len_max readLength(uint_reads_cnt_max i) { return reads[i].length(); };

            uint_read_len_max maxReadLengthVirtual() { return maxReadLength(); };
            uint_reads_cnt_max readsCountVirtual() { return readsCount(); };
            bool isReadLengthConstantVirtual() { return isReadLengthConstant(); };
            const string getReadVirtual(uint_reads_cnt_max i) { return getRead(i); };
            uint_read_len_max readLengthVirtual(uint_reads_cnt_max i) { return readLength(i); };

            void printout() {
                properties->printout();

//                for (auto readIt = reads.begin(); readIt != reads.end(); readIt++)
//                    cout << *readIt << "\n";
            }
           
            static DefaultReadsSet* readReadsSet(std::string filename, std::string pairfile = "") {
                istream* streamSource = new std::ifstream(filename, std::ios::in | std::ios::binary);
                istream* pairSource = 0;
                if (pairfile != "")
                    pairSource = new std::ifstream(pairfile, std::ios::in | std::ios::binary);
                        
                DefaultReadsSet* readsSet;
                if (filename.substr(filename.length() - 6) == ".fasta") {
                    FASTAReadsSourceIterator<uint_read_len_max>* readsSource = new FASTAReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new DefaultReadsSet(readsSource);
                } else if (filename.substr(filename.length() - 6) == ".fastq") {
                    FASTQReadsSourceIterator<uint_read_len_max>* readsSource = new FASTQReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new DefaultReadsSet(readsSource);
                } else {
                    ConcatenatedReadsSourceIterator<uint_read_len_max>* readsSource = new ConcatenatedReadsSourceIterator<uint_read_len_max>(streamSource);
                    readsSet = new DefaultReadsSet(readsSource);
                }
                delete(streamSource);
                return readsSet;
            }
    };

    //typedef DefaultReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_max>> ConcatenatedReadsSet;

    template<class ReadsSourceIterator>
    DefaultReadsSet::DefaultReadsSet(ReadsSourceIterator* readsIterator) {

        bool symbolOccured[UCHAR_MAX] = {0};

        while (readsIterator->moveNext()) {

            properties->readsCount++;

            // analize read length
            uint_read_len_max length = readsIterator->getReadLength();
            if (properties->maxReadLength == 0)
                properties->maxReadLength = length;
            else if (properties->maxReadLength != length) {
                properties->constantReadLength = false;
                if (properties->maxReadLength < length)
                    properties->maxReadLength = length;
            }

            properties->allReadsLength += length;

            //analize symbols
            string read(readsIterator->getRead());

            for (uint_read_len_max i = 0; i < length; i++) {
                read[i] = toupper(read[i]);
                if (!symbolOccured[(unsigned char) read[i]]) {
                    symbolOccured[(unsigned char) read[i]] = true;
                    properties->symbolsCount++;
                }
            }

            reads.push_back(read);
        }

        // order symbols

        for (uint_symbols_cnt i = 0, j = 0; i < properties->symbolsCount; j++)
            if (symbolOccured[(unsigned char) j])
                properties->symbolsList[(unsigned char) (i++)] = j;

        properties->generateSymbolOrder();

    };

}

#endif // DEFAULTREADSSET_H_INCLUDED
