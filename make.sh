#!/usr/bin/env bash

g++ -E findIndicatorMats.cpp -o main.i
g++ -S main.i -o main.s
g++ -c main.s -o main.o
g++ -std=c++23 main.o -o main
