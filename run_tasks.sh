#!/bin/bash

FILES=./task/*

main()
{
    for f in $FILES
    do
        #echo $f
        build/tsp $f
    done
}

main
