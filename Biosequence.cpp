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

size_t Biosequence::size() {
    return len;
}

size_t Biosequence::datastore_size() {
    return datastore.size();
}

Biosequence::Biosequence(unsigned int n) {
    datastore = std::vector<uint64_t>(seq_datastore_len(n));
    len = n;
}

SeqIndexer Biosequence::operator[](unsigned int idx) {
    return SeqIndexer(*this, idx);
}

int64_t SeqIndexer::index() {
    return position >> 6;
}

int64_t SeqIndexer::offset() {
    return position & 0b111111;
}