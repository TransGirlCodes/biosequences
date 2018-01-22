//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "Biosequence.h"
#include "bitops.h"

uint64_t seq_datastore_len(unsigned int len) {
    auto nucs_per_chunk = 64 / 4;
    auto d = len / nucs_per_chunk;
    return d + (((len > 0) == (nucs_per_chunk > 0)) & (d * nucs_per_chunk != len));
}

// Biosequence methods.

const Biosequence::size_type Biosequence::size() const {
    return len;
}

const Biosequence::size_type Biosequence::datastore_size() const {
    return datastore.size();
}

Biosequence::Biosequence(size_type n) {
    datastore = std::vector<uint64_t>(seq_datastore_len(n));
    len = n;
}

Biosequence::size_type Biosequence::reference::index() const {
    return position >> 6;
}

Biosequence::size_type Biosequence::reference::offset() const {
    return position & 0b111111;
}