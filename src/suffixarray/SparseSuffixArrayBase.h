/* 
 * File:   SparseSuffixArrayBase.h
 * Author: Tomek
 *
 * Created on 18 luty 2014, 10:01
 */

#ifndef SPARSESUFFIXARRAYBASE_H
#define	SPARSESUFFIXARRAYBASE_H

#include "SuffixArrayBase.h"

namespace PgSAIndex {

    class SparseSuffixArrayBase: public SuffixArrayBase
    {
        protected:
            const uchar skippedSymbolsCount; 
            const uchar bytesPerElement; 
            
            SparseSuffixArrayBase(uint_max elementsCount, uint_max suffixArrayBytes, PseudoGenomeBase* pg, uchar skippedSymbolsCount, uchar bytesPerElement)
            : SuffixArrayBase(elementsCount, suffixArrayBytes, pg),
              skippedSymbolsCount(skippedSymbolsCount),
              bytesPerElement(bytesPerElement)
            { };

            virtual ~SparseSuffixArrayBase() {};

        public:

            const uchar getSkippedSymbolsCount() { return this->skippedSymbolsCount; };
            const uchar getBytesPerElement() { return this->bytesPerElement; };
      
            bool isSAElementMinimal() { return bytesPerElement == 1; };
            bool isSAElementStandard() { return bytesPerElement == 2; };
    };


    class SparseSuffixArrayHeaderExtension {
        private:
            uchar skippedSymbolsCount;
            uchar bytesPerElement;

        public:

            SparseSuffixArrayHeaderExtension(SparseSuffixArrayBase* ssab)
            : skippedSymbolsCount(ssab->getSkippedSymbolsCount()),
              bytesPerElement(ssab->getBytesPerElement())
            {}

            SparseSuffixArrayHeaderExtension(std::istream& src) {
                int srchelper;
                src >> srchelper;
                skippedSymbolsCount = srchelper;
                src >> srchelper;
                bytesPerElement = srchelper;                
                src.get();
            }

            ~SparseSuffixArrayHeaderExtension() { };
            
            void write(std::ostream& dest) {
                dest << (int) skippedSymbolsCount << "\n";
                dest << (int) bytesPerElement << "\n";
            }

            bool isSAElementMinimal() { return bytesPerElement == 1; };
            bool isSAElementStandard() { return bytesPerElement == 2; };
    };
}

#endif	/* SPARSESUFFIXARRAYBASE_H */

