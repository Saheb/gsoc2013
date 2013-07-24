#!/bin/sh
g++ -I../include -L../_release -O0 -omain "$@" -lOGDF -lCOIN -pthread
