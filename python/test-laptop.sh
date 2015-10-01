#!/bin/bash

./clear.sh
./run.sh \
    -m 1024 \
    -s 2 \
    -r 1000 \
    -R True \
    -S $$ \
    -c 2 \
    2>&1

