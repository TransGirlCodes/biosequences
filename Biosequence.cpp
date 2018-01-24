//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "Biosequence.h"

uint64_t seq_datastore_len(unsigned long len, int bits_per_elem) {
    auto symbols_per_chunk = 64 / bits_per_elem;
    auto d = len / symbols_per_chunk;
    return d + (((len > 0) == (symbols_per_chunk > 0)) & (d * symbols_per_chunk != len));
}

// Concrete algorithm for the complement of a sequence encoded using the IUPAC alphabet.
template<>
void Biosequence<IUPAC_DNA>::complement() {
    auto next = reference(*this, 0);
    auto stop = reference(*this, len);
    while(next < stop) {
        auto chunk = next.get_chunk();
        auto new_chunk = (((chunk & 0x1111111111111111) << 3) | ((chunk & 0x8888888888888888) >> 3) |
                          ((chunk & 0x2222222222222222) << 1) | ((chunk & 0x4444444444444444) >> 1));
        next.insert_chunk(new_chunk);
        next.increment_by(64);
    }
}

// Concrete algorithm for the complement of a sequence encoded using the ATCG alphabet.
template<>
void Biosequence<ATCG>::complement() {
    auto next = reference(*this, 0);
    auto stop = reference(*this, len);
    while(next < stop) {
        auto chunk = next.get_chunk();
        auto new_chunk = ~chunk;
        next.insert_chunk(new_chunk);
        next.increment_by(64);
    }
}

template<typename A>
Biosequence<A> Biosequence<A>::reverse() {
    auto next = reference(*this, len);
    auto stop = reference(*this, 0);
    auto i = 1;
    auto r = next.bits_remaining();
}







