#!/bin/bash

for i in Bernoulli Uniform Normal Mixture; do
    for j in QLearning UCB1 ParticleFilter Beta NormalGamma DirichletNormalGamma; do
        rm -f $i-$j:*
    done
done

rm -f MixtureBandit-*
rm -f curve.png
rm -f core*
rm -f gnuplot_cmds
