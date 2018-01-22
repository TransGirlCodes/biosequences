//
// Created by Ben Ward on 25/09/2017.
//

#include <type_traits>
#include <numeric>
#include <iostream>
#include "../biosymbols.h"

int main(){

    biosymbols::DNA a = biosymbols::DNA::A;
    char b = biosymbols::as_character(biosymbols::DNA::T);
    std::cout << "The DNA symbol is: " << a << std::endl;
    std::cout << "The symbol is: " << b << std::endl;
    auto x = biosymbols::DNA('G');
    std::cout << biosymbols::as_character(~x) << std::endl;

    return 0;
}
