#!/bin/bash
mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8
pwd
make test
(touch ./docs/html/.nojekyll || true)