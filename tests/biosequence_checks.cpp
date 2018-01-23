//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "../biosequence.h"
#include "gtest/gtest.h"

// Biosymbols.
using namespace biosymbols;

TEST(Biosymbols, DNAIntoChar) {
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


TEST(SequenceConstruction, EmptySequences) {
    EXPECT_EQ(26, Biosequence(26).size());
}