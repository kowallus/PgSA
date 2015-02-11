#ifndef SUFFIXARRAYLOOKUPTABLEINTERFACE_H_INCLUDED
#define SUFFIXARRAYLOOKUPTABLEINTERFACE_H_INCLUDED

#include "../SuffixArrayBase.h"

namespace PgSAIndex {

    template <typename uint_read_len, typename uint_pg_len>
    class SuffixArrayLookupTableInterface
    {
        private:
        protected:
            uint_read_len const keyPrefixLength;

            SuffixArrayLookupTableInterface(uint_read_len keyPrefixLength)
            : keyPrefixLength(keyPrefixLength)
            {}

        public:
            virtual ~SuffixArrayLookupTableInterface() {};

            inline uint_read_len getKeyPrefixLength() { return keyPrefixLength; }
            
            virtual void write(std::ostream& dest) = 0;

    };
}

#endif // SUFFIXARRAYLOOKUPTABLEINTERFACE_H_INCLUDED
