#!/bin/sh

PROJECT="."
PROJECT_ROOT="$PROJECT"

SOURCE_DIR="$PROJECT_ROOT/src"
INCLUDE_DIR="$PROJECT_ROOT/include"

# compile
clang -c -O3 -Wall -fpic -march=native -std=c++17 -mfma $SOURCE_DIR/*.cpp \
  -I $INCLUDE_DIR


# # build library
g++ *.o -o test.out
g++ -shared *.o -o _nig.so

# move objects to build
mkdir -p build
mv *.o build/