//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include <vector>
#include "biosymbols.h"
#include "bitops.h"

#ifndef BIOSEQUENCES_LIBRARY_H
#define BIOSEQUENCES_LIBRARY_H

// Stand alone function, calculate how big a Biosequence's
// data store should be.
uint64_t seq_datastore_len(unsigned int len);

class Biosequence {

    typedef biosymbols::DNA value_type;
    typedef uint64_t size_type;

private:
    std::vector<uint64_t> datastore;
    size_type len;

public:

    // Individual symbol reference.
    class reference {

    private:
        size_type position = 0;
        Biosequence& sequence;
        uint64_t get_chunk() const {
            return sequence.datastore[index()];
        }

    public:
        reference(Biosequence& seq, size_type pos)
        : sequence(seq)
        {
            position = pos << trailing_zeros((size_type) 4);
        }
        size_type index() const;
        size_type offset() const;

        // Converts the reference into a value.
        operator value_type() const {
            // Extract chunk of data.
            uint64_t chunk = get_chunk();
            // Shift chunk of data, you now have the encoded symbol.
            uint64_t encoded_symbol = (chunk >> offset()) & uint64_t(15);
            // TODO: Add alphabet dependent decoder step here.
            // NOTE: With 4 bit DNA alphabet there isn't any magic
            // needed to decode the symbol from the binary.
            return static_cast<biosymbols::DNA>(encoded_symbol);
        }

        // Set the value of a nucleotide in a DNA sequence.
        reference& operator=(value_type x) {
            // Convert the symbol to a uint64_t.
            auto encoded_symbol = static_cast<uint64_t>(x);
            uint64_t offset_symbol = encoded_symbol << offset();

            uint64_t chunk_mask = ~(uint64_t(15) << offset());
            uint64_t masked_chunk = get_chunk() & chunk_mask;
            uint64_t new_chunk = offset_symbol | masked_chunk;
            sequence.datastore[index()] = new_chunk;
            return *this;
        }

        value_type operator~() const {
            value_type element = *this;
            return ~element;
        }

    };

    Biosequence(size_type n);
    const size_type size() const;
    const size_type datastore_size() const;
    reference operator[]( size_type pos ) {
        return reference(*this, pos);
    };

};

std::ostream& operator<<(std::ostream &out, Biosequence& seq) {
    for(int i = 0; i < seq.size(); i++) {
        out << biosymbols::as_character(biosymbols::DNA(seq[i]));
    }
    return out;
}

#endif