#!/bin/bash

set -e

TEXFILE=rocket.tex

CONF=$1
if [ -z $CONF ]; then
    echo "no config given"
    exit 1
fi

./bin/prog/rocketgen -c $CONF > $TEXFILE

pdflatex $TEXFILE
