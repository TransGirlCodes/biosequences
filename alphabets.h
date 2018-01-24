//
// Created by Ben Ward (EI) on 22/01/2018.
//

#include <exception>

#ifndef BIOSEQUENCES_ALPHABETS_H
#define BIOSEQUENCES_ALPHABETS_H

class IUPAC_DNA {
    typedef biosymbols::DNA symbol_type;
    static const std::string alphabet_name;
public:
    static const int bits_per_symbol = 4;

    static uint64_t encode(symbol_type nt) {
        return static_cast<uint64_t>(nt);
    }

    static symbol_type decode(uint64_t binary) {
        return static_cast<symbol_type>(binary);
    }
};

const std::string IUPAC_DNA::alphabet_name = "IUPAC DNA Alphabet";


#endif //BIOSEQUENCES_ALPHABETS_H
