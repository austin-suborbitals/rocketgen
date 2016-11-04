#!/bin/bash

set -e

OPTS=-d
CONVTOOL=rubber

CONF=$1
if [ -z $CONF ]; then
    echo "no config given"
    exit 1
fi

OUTNAME=rocket
if [ "$#" -eq 2 ]; then
    OUTNAME=$2
fi

TEXFILE=$OUTNAME.tex

./bin/prog/rocketgen -c $CONF > $TEXFILE

$CONVTOOL $OPTS $TEXFILE
