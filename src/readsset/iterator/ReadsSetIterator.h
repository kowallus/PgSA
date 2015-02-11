#ifndef READSSETITERATOR_H_INCLUDED
#define READSSETITERATOR_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include "../../helper.h"
#include "../../pgsaconfig.h"

namespace PgSAReadsSet {

    typedef uint_read_len_std uint_read_len_max;

    template < typename uint_read_len >
    class ReadsSourceIteratorTemplate
    {
        public:

            virtual ~ReadsSourceIteratorTemplate() {};

            virtual bool moveNextVirtual() = 0;
            virtual string getReadVirtual() = 0;
            virtual uint_read_len getReadLengthVirtual() = 0;
    };

    template < typename uint_read_len >
    class ConcatenatedReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;

        public:

            ConcatenatedReadsSourceIterator(std::istream* source) {
                this->source = source;
            }

            ~ConcatenatedReadsSourceIterator() {};

            bool moveNext() {
                if (!std::getline(*source, line))
                    return false;

                for(length = 0; length < line.length(); length++)
                    if (!isalpha(line[length]))
                        break;

                return true;
            };

            string getRead() {
                return line.substr(0, length);
            };

            uint_read_len getReadLength() {
                return length;
            };

            bool moveNextVirtual() { return moveNext(); };
            string getReadVirtual() { return getRead(); };
            uint_read_len getReadLengthVirtual() { return getReadLength(); };

    };
    
    template < typename uint_read_len >
    class FASTAReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;
            std::istream* pairSource = 0;
            bool pair = false;
            
        public:

            FASTAReadsSourceIterator(std::istream* source) {
                this->source = source;
            }
            
            FASTAReadsSourceIterator(std::istream* source, std::istream* pairSource) {
                this->source = source;
                this->pairSource = pairSource;
            }

            ~FASTAReadsSourceIterator() {};

            bool moveNext() {
                std::istream* src = source;
                if (pair && pairSource)
                    src = pairSource;
                pair = !pair;
                do {
                    if (!std::getline(*src, line))
                        return false;
                } while (line.find('>') == 0);

                for(length = 0; length < line.length(); length++)
                    if (!isalpha(line[length]))
                        break;

                return true;
            };

            string getRead() {
                return line.substr(0, length);
            };

            uint_read_len getReadLength() {
                return length;
            };

            bool moveNextVirtual() { return moveNext(); };
            string getReadVirtual() { return getRead(); };
            uint_read_len getReadLengthVirtual() { return getReadLength(); };
    };
    
    template < typename uint_read_len >
    class FASTQReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;
            std::istream* pairSource = 0;
            bool pair = false;
            
        public:
            
            FASTQReadsSourceIterator(std::istream* source, std::istream* pairSource = std::string()) {
                this->source = source;
                this->pairSource = pairSource;
            }

            ~FASTQReadsSourceIterator() {};

            bool moveNext() {
                std::istream* src = source;
                if (pair && pairSource)
                    src = pairSource;
                pair = !pair;
                
                string someinfo;
                if (!std::getline(*src, someinfo))
                    return false;
                std::getline(*src, line);
                std::getline(*src, someinfo);
                std::getline(*src, someinfo);

                for(length = 0; length < line.length(); length++)
                    if (!isalpha(line[length]))
                        break;

                return true;
            };

            string getRead() {
                return line.substr(0, length);
            };

            uint_read_len getReadLength() {
                return length;
            };

            bool moveNextVirtual() { return moveNext(); };
            string getReadVirtual() { return getRead(); };
            uint_read_len getReadLengthVirtual() { return getReadLength(); };
    };
}

#endif // READSSETITERATOR_H_INCLUDED
