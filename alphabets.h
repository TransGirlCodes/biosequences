//
// Created by Ben Ward (EI) on 22/01/2018.
//

#include <exception>

#ifndef BIOSEQUENCES_ALPHABETS_H
#define BIOSEQUENCES_ALPHABETS_H

class IUPAC_DNA {
    static const std::string alphabet_name;
public:
    typedef biosymbols::DNA symbol_type;
    static const int bits_per_symbol = 4;

    static uint64_t encode(symbol_type nt) {
        return static_cast<uint64_t>(nt);
    }

    static symbol_type decode(uint64_t binary) {
        return static_cast<symbol_type>(binary);
    }
};

const std::string IUPAC_DNA::alphabet_name = "IUPAC DNA Alphabet";

class ATCG {
    static const std::string alphabet_name;
public:
    typedef biosymbols::DNA symbol_type;
    static const int bits_per_symbol = 2;

    static uint64_t encode(symbol_type nt) {
        return (uint64_t) trailing_zeros(nt);

    };

    static symbol_type decode(uint64_t binary) {
        return static_cast<symbol_type>(uint64_t(1) << binary);
    };
};




#endif //BIOSEQUENCES_ALPHABETS_H
