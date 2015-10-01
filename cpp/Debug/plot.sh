#!/bin/bash

GNUPLOT_SCRIPT="gnuplot_cmds"

gnuplot $GNUPLOT_SCRIPT
eog curve.png
