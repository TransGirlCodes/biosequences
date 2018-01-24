//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include <vector>
#include <string>
#include "biosymbols.h"
#include "bitops.h"
#include "alphabets.h"

#ifndef BIOSEQUENCES_LIBRARY_H
#define BIOSEQUENCES_LIBRARY_H

// Stand alone function, calculate how big a Biosequence's
// data store should be.
uint64_t seq_datastore_len(unsigned long len);


/*! \class Biosequence
 *  \brief A representation of a biological sequence.
 *
 *  The Biosequence class compactly stores a sequence, which can be
 *  thought of as a vector of DNA or RNA symbols.
 *
 *  Symbols are stored in the datastore of a Biosequence using some
 *  encoding that is defined by an alphabet class.
 *
 * For a given nucleotide position j, Biosequence::reference indexes into
 * data store to extract the nucleotide.
 *
 * index :        index(j) - 1       index(j)       index(j) + 1
 * data : ....|xxxxx...........|xxXxxxxxxxxxxxxx|............xxxx|....
 * offset :                        |<-offset(j)-|
 * width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
 *
 * '.' : unused (4 bits/char)
 * 'x' : used
 * 'X' : used and pointed by index `i`
 */
class Biosequence {

    typedef typename biosymbols::DNA value_type;
    typedef uint64_t size_type;

private:
    //! A vector holding the binary encoded sequence representation.
    std::vector<uint64_t> datastore;
    //! How long the sequence is, in number of nucleotides.
    size_type len = 0;

public:

    // Individual symbol reference.
    class reference {
    private:
        size_type position = 0;
        Biosequence& sequence;

        uint64_t get_chunk() const {
            return sequence.datastore[index()];
        }

        uint64_t extract_encoded_symbol() const {
            auto chunk = get_chunk();
            return (chunk >> offset()) & uint64_t(15);
        }

    public:
        reference(Biosequence& seq, size_type pos)
        : sequence(seq)
        {
            position = pos << bitops::trailing_zeros((size_type) 4);
        }
        size_type index() const;
        size_type offset() const;

        // Converts the reference into a value.
        operator value_type() const {
            auto encoded_symbol = extract_encoded_symbol();
            return IUPAC_DNA::decode(encoded_symbol);
        }

        // Set the value of a nucleotide in a DNA sequence, with a symbol.
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

        reference& operator++() { // prefix increment.
            position += IUPAC_DNA::bits_per_symbol;
            return *this;
        }

        reference operator++ (int) { // postfix increment.
            reference result(*this);
            ++(*this);
            return result;
        }


    };

    Biosequence(size_type n);

    Biosequence(std::string& str) {
        len = str.size();
        datastore.resize(seq_datastore_len(str.size()));
        reference i(*this, 0);
        for(char& c : str) {
            i = biosymbols::as_symbol<biosymbols::DNA>(c);
            i++;
        }
    }

    const size_type size() const;

    const size_type datastore_size() const;

    reference operator[](size_type pos) {
        return reference(*this, pos);
    };

    operator std::string() {
        auto str = std::string(size(), 'X');
        auto r = reference(*this, 0);
        for(int i = 0; i < size(); i++) {
            char letter = biosymbols::as_character((value_type) r);
            str[i] = letter;
            r++;
        }
        return str;
    }

};

std::ostream& operator<<(std::ostream &out, Biosequence& seq) {
    for(int i = 0; i < seq.size(); i++) {
        out << biosymbols::as_character((biosymbols::DNA) seq[i]);
    }
    return out;
}

#endif