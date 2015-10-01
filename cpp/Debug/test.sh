#! /bin/bash

SEED=$RANDOM
GNUPLOT_SCRIPT="gnuplot_cmds"
BANDIT="Mixture"
BASIS="4"
ARMS2="3"
TRIALS2="10"
RUNS2="10"

check_id() {
    cat ../src/utils.h | grep "\<$1\s=" | awk '{print $3}' | sed "s/,//g"
}

run() {
    AGENT="$1"
    AGENT_ID=`check_id $AGENT`
    
    shift
    OUTPUT="$BANDIT-$AGENT:$*"

    echo -n " '$OUTPUT' w l," >> $GNUPLOT_SCRIPT
    echo "./bandit --seed $SEED --bandit $BANDIT_ID --basis $BASIS --agent $AGENT_ID --arms $ARMS2 --trials $TRIALS2 --runs $RUNS2 $* > $OUTPUT"
    time ./bandit --seed $SEED --bandit $BANDIT_ID --basis $BASIS --agent $AGENT_ID --arms $ARMS2 --trials $TRIALS2 --runs $RUNS2 $* > "$OUTPUT" &
}

if [ ! -z $1 ]; then
    SEED=$1
fi

BANDIT_ID=`check_id $BANDIT`

./clear.sh

cat << EOF > $GNUPLOT_SCRIPT
set terminal png font '/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansCondensed.ttf,11' size 1280,800
set output "curve.png"
set grid
set log x
EOF

echo -n "plot " >> $GNUPLOT_SCRIPT

run UCB1 --exploration 1
run NormalGamma --beta 1
run DirichletNormalGamma --beta 1
run Beta
run QLearning --eplison 0.1

wait

sed -i -e 's/,$//g' $GNUPLOT_SCRIPT
./plot.sh &
 
