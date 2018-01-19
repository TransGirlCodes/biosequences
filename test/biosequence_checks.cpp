//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "../biosequence.h"
#include <iostream>

int main(){

    Biosequence seq = Biosequence(26);

    std::cout << "Sequence has size: " << seq.size() << std::endl;
    std::cout << "Sequence has datasize: " << seq.datastore_size() << std::endl;

    SeqIndexer a = seq[12];
    SeqIndexer b = seq[15];
    SeqIndexer c = seq[16];

    std::cout << "a has index of: " << a.index() << "\nand an offset of: " <<
            a.offset() << std::endl;

    std::cout << "b has index of: " << b.index() << "\nand an offset of: " <<
              b.offset() << std::endl;

    std::cout << "c has index of: " << c.index() << "\nand an offset of: " <<
              c.offset() << std::endl;
    return 0;
}