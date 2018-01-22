//
// Created by Ben Ward (EI) on 18/01/2018.
//

#include "../biosequence.h"
#include "gtest/gtest.h"

TEST(SequenceConstruction, EmptySequences) {
    EXPECT_EQ(26, Biosequence(26).size());
}