#!/bin/bash

input=$1

softpoint.x leshouches < $input > slhaspectrum.in
./run

