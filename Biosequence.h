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
uint64_t seq_datastore_len(unsigned long len, int bits_per_elem);

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
template <typename A>
class Biosequence {

    typedef typename A::symbol_type value_type;
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

        uint64_t symbol_mask() const {
            return (uint64_t(1) << A::bits_per_symbol) - uint64_t(1);
        }

        uint64_t get_masked_chunk() const {
            uint64_t chunk_mask = ~(symbol_mask() << offset());
            return get_chunk() & chunk_mask;
        }

        uint64_t extract_encoded_symbol() const {
            auto chunk = get_chunk();
            return (chunk >> offset()) & uint64_t(15);
        }

        void insert_encoded_symbol(uint64_t encoded_symbol) {
            uint64_t offset_symbol = encoded_symbol << offset();
            auto old_masked_chunk = get_masked_chunk();
            uint64_t new_chunk = offset_symbol | old_masked_chunk;
            insert_chunk(new_chunk);
        }

    public:
        reference(Biosequence& seq, size_type pos)
        : sequence(seq)
        {
            position = pos << bitops::trailing_zeros((size_type) 4);
        }

        size_type index() const {
            return position >> 6;
        };

        size_type offset() const {
            return position & 0b111111;
        };

        uint64_t get_chunk() const {
            return sequence.datastore[index()];
        }

        void insert_chunk(uint64_t chunk) {
            sequence.datastore[index()] = chunk;
        }

        // Converts the reference into a value.
        operator value_type() const {
            auto encoded_symbol = extract_encoded_symbol();
            return A::decode(encoded_symbol);
        }

        // Set the value of a nucleotide in a DNA sequence, with a symbol.
        reference& operator=(value_type x) {
            auto encoded_symbol = A::encode(x);
            insert_encoded_symbol(encoded_symbol);
            return *this;
        }

        value_type operator~() const {
            value_type element = *this;
            return ~element;
        }

        reference& operator++() { // prefix increment.
            position += A::bits_per_symbol;
            return *this;
        }

        reference operator++(int) { // postfix increment.
            reference result(*this);
            ++(*this);
            return result;
        }

        reference& operator--() {  // prefix decrement.
            position -= A::bits_per_symbol;
            return *this;
        }

        reference operator--(int) { // postfix increment.
            reference result(*this);
            --(*this);
            return result;
        }

        int bits_remaining() {
            return (offset() + A::bits_per_symbol) % 64;
        }

        void increment_by(int x) {
            position += x;
        }

        friend bool operator==(const reference& r1, const reference& r2) {
            return (r1.position == r2.position) && (&r1.sequence == &r2.sequence);
        }

        friend bool operator<(const reference& r1, const reference& r2) {
            // TODO: Throw an exception if the two references do no point to a biological sequence.
            return r1.position < r2.position;
        }

    };

    explicit Biosequence(size_type n) {
        datastore = std::vector<uint64_t>(seq_datastore_len(n, A::bits_per_symbol));
        len = n;
    }

    explicit Biosequence(std::string& str) {
        len = str.size();
        datastore.resize(seq_datastore_len(str.size(), A::bits_per_symbol));
        reference i(*this, 0);
        for(char& c : str) {
            i = biosymbols::as_symbol<value_type>(c);
            i++;
        }
    }

    const size_type size() const {
        return len;
    };

    const size_type datastore_size() const {
        return datastore.size();
    };

    reference operator[](size_type pos) {
        return reference(*this, pos);
    };

    explicit operator std::string() {
        auto str = std::string(size(), 'X');
        auto r = reference(*this, 0);
        for(int i = 0; i < size(); i++) {
            char letter = biosymbols::as_character((value_type) r);
            str[i] = letter;
            r++;
        }
        return str;
    }

    friend std::ostream& operator<<(std::ostream &output, Biosequence<A>& seq) {
        for(int i = 0; i < seq.size(); i++) {
            output << biosymbols::as_character((value_type) seq[i]);
        }
        return output;
    }

    void complement();

    Biosequence<A> reverse();

};

#endif