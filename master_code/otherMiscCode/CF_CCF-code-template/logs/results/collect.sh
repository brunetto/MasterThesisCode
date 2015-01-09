#!/usr/bin/env bash

cat ../*-result.dat | grep Particles > particles.dat
cat ../*-result.dat | grep Random > random.dat
cat ../*-result.dat | grep Mixed_counts_1 > m1.dat
cat ../*-result.dat | grep Mixed_counts_2 > m2.dat
cat ../*-result.dat | grep Radii > radii.dat