#!/bin/bash

FILES=./task/*

main()
{
    for f in $FILES
    do
        build/tsp $f
    done
}

main
