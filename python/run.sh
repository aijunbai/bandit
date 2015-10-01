#!/bin/bash

PYTHON=`which pypy`

INP="bandit.txt"
MIN="1"
MAX="1024"
STP="2"
REP="1000"
RND="True"
SED="$$"

CPU=`cat /proc/cpuinfo | grep processor | wc -l`
NUM=(8 16 32 64 128 256 512)
ALG=(RoundRobin Randomized Greedy UCB ThompsonSampling)

while getopts  "m:s:r:R:S:c:" flag; do
    case "$flag" in
        m) MAX=$OPTARG;;
        s) STP=$OPTARG;;
        r) REP=$OPTARG;;
        R) RND=$OPTARG;;
        S) SED=$OPTARG;;
        c) CPU=$OPTARG;;
    esac
done

run() {
    local NUM=$1
    local ALG=$2
    local DIR="output-${NUM}x${REP}.dir"
    local OUT="output-$ALG.txt"

    mkdir -p $DIR
    cd $DIR
    rm -f $OUT

    $PYTHON ../experiment.py \
        --min $MIN --max $MAX --step $STP --repeat $REP --random $RND --bandits $NUM --seed $SED $INP $ALG \
        2>&1 1>$OUT

    cd ..
}

job() {
    local ID="$1"

    for i in `seq 1 ${#NUM[@]}`; do
        for j in `seq 1 ${#ALG[@]}`; do
            local num=`expr $i - 1`
            local alg=`expr $j - 1`
            local task=`expr $num \* ${#ALG[@]} + $j`

            if [ `expr $task % $CPU` -eq $ID ]; then
                echo running ${NUM[$num]} ${ALG[$alg]} @ $ID...
                run ${NUM[$num]} ${ALG[$alg]}
                echo finished ${NUM[$num]} ${ALG[$alg]} @ $ID...
            fi
        done
    done
}

for i in `seq 1 $CPU`; do
    job `expr $i - 1` &
done

