#!/bin/bash

./clear.sh
./run.sh \
    -m 32768 \
    -s 2 \
    -r 10000 \
    -R True \
    -S $$ \
    -c 7 \
    2>&1

