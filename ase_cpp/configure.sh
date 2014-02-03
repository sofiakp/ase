#!/bin/bash -x

mkdir -p ./build/release
mkdir -p ./build/debug
mkdir -p ./bin/debug

cd build/release
cmake -D CMAKE_BUILD_TYPE:STRING=release ../..
cd ../..

cd build/debug
cmake -D CMAKE_BUILD_TYPE:STRING=debug ../..
cd ../..

