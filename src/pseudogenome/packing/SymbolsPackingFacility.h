/* 
 * File:   SymbolsPackingFacility.h
 * Author: coach
 *
 * Created on 12 luty 2014, 11:43
 */

#ifndef SYMBOLSPACKINGFACILITY_H
#define	SYMBOLSPACKINGFACILITY_H

#include "../../helper.h"
#include "../../pgsaconfig.h"
#include "../../readsset/ReadsSetBase.h"

using namespace PgSAHelpers;

namespace PgSAIndex {
    
    template<typename uint_element>
    class SymbolsPackingFacility {
        private:
            
            uint_max maxValue;
            const uint_symbols_cnt symbolsCount;
            const uchar symbolsPerElement;
            
            char symbolsList[UCHAR_MAX] = {0};
            int symbolOrder[UCHAR_MAX] = {-1};

            // for a given value reverse[value][pos] returns character at the given position
            char_pg** reverse;
            char_pg* reverseFlat;
            
            // for a given value clear[value][prefixLenght] returns value of first prefixLength symbols followed by 0 symbols.
            uint_element** clear;
            uint_element* clearFlat;
            
            void buildReverseAndClearIndexes();
            
        public:
            SymbolsPackingFacility(ReadsSetProperties* readsSetProperties, uchar symbolsPerElement);
            
            ~SymbolsPackingFacility();
            
            uint_element getMaxValue();
            
            // value should be not greater then maxValue
            char_pg reverseValue(uint_element value, uchar position);
            
            // value should be not greater then maxValue
            string reverseValue(uint_element value);
            
            const string reverseSequence(uint_element* sequence, const uint_max pos, const uint_max length);
            
            void reverseSequence(uint_element* sequence, const uint_max pos, const uint_max length, char_pg* kmerPtr);  
          
            // sequence should consist of at least symbolsPerElement symbols; result should be not greater then maxValue
            uint_element packSymbols(const char_pg* symbols);
            
            // result should be not greater then maxValue
            uint_element packSuffixSymbols(const char_pg* symbols, const uint_max length);
            
            // result should be not greater then maxValue
            uint_element packPrefixSymbols(const char_pg* symbols, const uint_max length);
            
            uint_max packSequence(const char_pg* source, const uint_max length, uint_element* dest);

            uint_element clearSuffix(const uint_element value, uchar prefixLength);
            
            static bool isCompatibile(uchar symbolsPerElement, uchar symbolsCount);
            
            static uchar maxSymbolsPerElement(uchar symbolsCount);
          
    };
    
}

#endif	/* SYMBOLSPACKINGFACILITY_H */

