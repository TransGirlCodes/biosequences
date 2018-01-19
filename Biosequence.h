//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include <vector>
#include "biosymbols.h"
#include "bitops.h"

#ifndef BIOSEQUENCES_LIBRARY_H
#define BIOSEQUENCES_LIBRARY_H

// Forward declare the SeqIndexer iterator.
class SeqIndexer;

// Stand alone function, calculate how big a Biosequence's
// data store should be.
uint64_t seq_datastore_len(unsigned int len);

class Biosequence {
private:
    std::vector<uint64_t> datastore;
    unsigned int len;

public:
    Biosequence(unsigned int n);
    size_t size();
    size_t datastore_size();

    SeqIndexer operator[](unsigned int idx);
};


// SeqIndexer is a class to refer to a position inside a sequence.
// TODO: Make this safe!
class SeqIndexer {
private:
    Biosequence& sequence;
    int64_t position;
public:
    SeqIndexer(Biosequence& seq, int64_t idx)
            : sequence(seq)
    {
        position = idx << trailing_zeros((uint64_t) 4);
    }
    int64_t index();
    int64_t offset();
};

#endif