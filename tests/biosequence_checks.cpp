//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "../biosequence.h"
#include "gtest/gtest.h"
#include <string>

// Biosymbols.
using namespace biosymbols;

TEST(biosymbols, DNAIntoChar) {
    EXPECT_EQ(as_character(DNA::A), 'A');
    EXPECT_EQ(as_character(DNA::T), 'T');
    EXPECT_EQ(as_character(DNA::C), 'C');
    EXPECT_EQ(as_character(DNA::G), 'G');
    EXPECT_EQ(as_character(DNA::Gap), '-');
    EXPECT_EQ(as_character(DNA::M), 'M');
    EXPECT_EQ(as_character(DNA::G), 'G');
    EXPECT_EQ(as_character(DNA::R), 'R');
    EXPECT_EQ(as_character(DNA::S), 'S');
    EXPECT_EQ(as_character(DNA::V), 'V');
    EXPECT_EQ(as_character(DNA::W), 'W');
    EXPECT_EQ(as_character(DNA::Y), 'Y');
    EXPECT_EQ(as_character(DNA::H), 'H');
    EXPECT_EQ(as_character(DNA::K), 'K');
    EXPECT_EQ(as_character(DNA::D), 'D');
    EXPECT_EQ(as_character(DNA::B), 'B');
    EXPECT_EQ(as_character(DNA::N), 'N');
}

TEST(biosymbols, CharIntoDNA) {
    EXPECT_EQ(as_symbol<DNA>('A'), DNA::A);
    EXPECT_EQ(as_symbol<DNA>('T'), DNA::T);
    EXPECT_EQ(as_symbol<DNA>('C'), DNA::C);
    EXPECT_EQ(as_symbol<DNA>('G'), DNA::G);
    EXPECT_EQ(as_symbol<DNA>('-'), DNA::Gap);
    EXPECT_EQ(as_symbol<DNA>('M'), DNA::M);
    EXPECT_EQ(as_symbol<DNA>('G'), DNA::G);
    EXPECT_EQ(as_symbol<DNA>('R'), DNA::R);
    EXPECT_EQ(as_symbol<DNA>('S'), DNA::S);
    EXPECT_EQ(as_symbol<DNA>('V'), DNA::V);
    EXPECT_EQ(as_symbol<DNA>('W'), DNA::W);
    EXPECT_EQ(as_symbol<DNA>('Y'), DNA::Y);
    EXPECT_EQ(as_symbol<DNA>('H'), DNA::H);
    EXPECT_EQ(as_symbol<DNA>('K'), DNA::K);
    EXPECT_EQ(as_symbol<DNA>('D'), DNA::D);
    EXPECT_EQ(as_symbol<DNA>('B'), DNA::B);
    EXPECT_EQ(as_symbol<DNA>('N'), DNA::N);
}

TEST(biosymbols, operations) {
    EXPECT_TRUE(is_gap(DNA::Gap));
    EXPECT_FALSE(is_gap(DNA::A));
    EXPECT_FALSE(is_gap(DNA::C));
    EXPECT_FALSE(is_gap(DNA::M));
    EXPECT_FALSE(is_gap(DNA::G));
    EXPECT_FALSE(is_gap(DNA::R));
    EXPECT_FALSE(is_gap(DNA::S));
    EXPECT_FALSE(is_gap(DNA::V));
    EXPECT_FALSE(is_gap(DNA::T));
    EXPECT_FALSE(is_gap(DNA::W));
    EXPECT_FALSE(is_gap(DNA::Y));
    EXPECT_FALSE(is_gap(DNA::H));
    EXPECT_FALSE(is_gap(DNA::K));
    EXPECT_FALSE(is_gap(DNA::D));
    EXPECT_FALSE(is_gap(DNA::B));
    EXPECT_FALSE(is_gap(DNA::N));

    EXPECT_FALSE(is_GC(DNA::Gap));
    EXPECT_FALSE(is_GC(DNA::A));
    EXPECT_TRUE(is_GC(DNA::C));
    EXPECT_FALSE(is_GC(DNA::M));
    EXPECT_TRUE(is_GC(DNA::G));
    EXPECT_FALSE(is_GC(DNA::R));
    EXPECT_TRUE(is_GC(DNA::S));
    EXPECT_FALSE(is_GC(DNA::V));
    EXPECT_FALSE(is_GC(DNA::T));
    EXPECT_FALSE(is_GC(DNA::W));
    EXPECT_FALSE(is_GC(DNA::Y));
    EXPECT_FALSE(is_GC(DNA::H));
    EXPECT_FALSE(is_GC(DNA::K));
    EXPECT_FALSE(is_GC(DNA::D));
    EXPECT_FALSE(is_GC(DNA::B));
    EXPECT_FALSE(is_GC(DNA::N));

    EXPECT_FALSE(is_purine(DNA::Gap));
    EXPECT_TRUE(is_purine(DNA::A));
    EXPECT_FALSE(is_purine(DNA::C));
    EXPECT_FALSE(is_purine(DNA::M));
    EXPECT_TRUE(is_purine(DNA::G));
    EXPECT_TRUE(is_purine(DNA::R));
    EXPECT_FALSE(is_purine(DNA::S));
    EXPECT_FALSE(is_purine(DNA::V));
    EXPECT_FALSE(is_purine(DNA::T));
    EXPECT_FALSE(is_purine(DNA::W));
    EXPECT_FALSE(is_purine(DNA::Y));
    EXPECT_FALSE(is_purine(DNA::H));
    EXPECT_FALSE(is_purine(DNA::K));
    EXPECT_FALSE(is_purine(DNA::D));
    EXPECT_FALSE(is_purine(DNA::B));
    EXPECT_FALSE(is_purine(DNA::N));

    EXPECT_FALSE(is_pyrimidine(DNA::Gap));
    EXPECT_FALSE(is_pyrimidine(DNA::A));
    EXPECT_TRUE(is_pyrimidine(DNA::C));
    EXPECT_FALSE(is_pyrimidine(DNA::M));
    EXPECT_FALSE(is_pyrimidine(DNA::G));
    EXPECT_FALSE(is_pyrimidine(DNA::R));
    EXPECT_FALSE(is_pyrimidine(DNA::S));
    EXPECT_FALSE(is_pyrimidine(DNA::V));
    EXPECT_TRUE(is_pyrimidine(DNA::T));
    EXPECT_FALSE(is_pyrimidine(DNA::W));
    EXPECT_TRUE(is_pyrimidine(DNA::Y));
    EXPECT_FALSE(is_pyrimidine(DNA::H));
    EXPECT_FALSE(is_pyrimidine(DNA::K));
    EXPECT_FALSE(is_pyrimidine(DNA::D));
    EXPECT_FALSE(is_pyrimidine(DNA::B));
    EXPECT_FALSE(is_pyrimidine(DNA::N));

    EXPECT_FALSE(is_ambiguous(DNA::Gap));
    EXPECT_FALSE(is_ambiguous(DNA::A));
    EXPECT_FALSE(is_ambiguous(DNA::C));
    EXPECT_TRUE(is_ambiguous(DNA::M));
    EXPECT_FALSE(is_ambiguous(DNA::G));
    EXPECT_TRUE(is_ambiguous(DNA::R));
    EXPECT_TRUE(is_ambiguous(DNA::S));
    EXPECT_TRUE(is_ambiguous(DNA::V));
    EXPECT_FALSE(is_ambiguous(DNA::T));
    EXPECT_TRUE(is_ambiguous(DNA::W));
    EXPECT_TRUE(is_ambiguous(DNA::Y));
    EXPECT_TRUE(is_ambiguous(DNA::H));
    EXPECT_TRUE(is_ambiguous(DNA::K));
    EXPECT_TRUE(is_ambiguous(DNA::D));
    EXPECT_TRUE(is_ambiguous(DNA::B));
    EXPECT_TRUE(is_ambiguous(DNA::N));

    EXPECT_TRUE(is_valid(DNA::Gap));
    EXPECT_TRUE(is_valid(DNA::A));
    EXPECT_TRUE(is_valid(DNA::C));
    EXPECT_TRUE(is_valid(DNA::M));
    EXPECT_TRUE(is_valid(DNA::G));
    EXPECT_TRUE(is_valid(DNA::R));
    EXPECT_TRUE(is_valid(DNA::S));
    EXPECT_TRUE(is_valid(DNA::V));
    EXPECT_TRUE(is_valid(DNA::T));
    EXPECT_TRUE(is_valid(DNA::W));
    EXPECT_TRUE(is_valid(DNA::Y));
    EXPECT_TRUE(is_valid(DNA::H));
    EXPECT_TRUE(is_valid(DNA::K));
    EXPECT_TRUE(is_valid(DNA::D));
    EXPECT_TRUE(is_valid(DNA::B));
    EXPECT_TRUE(is_valid(DNA::N));
    EXPECT_FALSE(is_valid(DNA::Invalid));
}

TEST(biosymbols, complement) {
    EXPECT_EQ(complement(DNA::A), DNA::T);
    EXPECT_EQ(complement(DNA::C), DNA::G);
    EXPECT_EQ(complement(DNA::G), DNA::C);
    EXPECT_EQ(complement(DNA::T), DNA::A);
    EXPECT_EQ(complement(DNA::Gap), DNA::Gap);
    EXPECT_EQ(complement(DNA::N), DNA::N);
}

TEST(biosymbols, arithmetic) {
    EXPECT_EQ(~DNA::Gap, DNA::N);
    EXPECT_EQ(DNA::Gap, ~DNA::N);
    EXPECT_EQ(trailing_zeros(DNA::A), 0);
    EXPECT_EQ(trailing_zeros(DNA::T), 3);
    EXPECT_EQ(trailing_zeros(DNA::C), 1);
    EXPECT_EQ(trailing_zeros(DNA::G), 2);
    EXPECT_EQ(trailing_zeros(DNA::N), 0);
    EXPECT_EQ(trailing_zeros(DNA::Gap), 8);
}

TEST(SequenceConstruction, EmptySequences) {
    EXPECT_EQ(26, Biosequence<IUPAC_DNA>(26).size());
}

TEST(biosequence, strings_and_printing) {
    // Test a round trip for making a seq from a string, and then a string from a sequence.
    std::ostringstream output;
    auto str = std::string("ATCGCCGAACGCGAAACGAAACCTGTGT");
    Biosequence<IUPAC_DNA> seq(str);
    std::string newstring(seq); // CLion flags this as error but it works fine with GCC.
    EXPECT_EQ(str, newstring);
    output << seq;
    EXPECT_EQ(output.str(), "ATCGCCGAACGCGAAACGAAACCTGTGT");

    auto complement_string = std::string("TAGCGGCTTGCGCTTTGCTTTGGACACA");
    seq.complement();
    std::string compstring(seq);
    EXPECT_EQ(complement_string, compstring);
}
