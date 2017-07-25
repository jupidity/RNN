#!/bin/bash

# shell script for building and cleaning nerual network
# Author: Sean Cassero


# compile with the proper directive
g++ -std=c++11 RNN.cpp -o main

# run the executable
./main

# clean the current directory
rm main
