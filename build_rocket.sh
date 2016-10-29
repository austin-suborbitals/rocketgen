#!/bin/bash

set -e

OPTS=-d
CONVTOOL=rubber

TEXFILE=rocket.tex

CONF=$1
if [ -z $CONF ]; then
    echo "no config given"
    exit 1
fi

./bin/prog/rocketgen -c $CONF > $TEXFILE

$CONVTOOL $OPTS $TEXFILE
