#!/bin/bash
# Create a synthetic burst from an image
# CFA images are built using raw = 1
# Standard burst: base_%i.tiff
# CFA images: base_cfa_%i.tiff
# Homographies: base_%i.hom

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH

# read input parameters
if [ "$#" -lt "5" ]; then
    echo "usage: $0 in base number sigma raw [interp boundary L type zoom]"
    echo "example: $0 im.tiff burst 100 5 0"
    exit 1
fi

IN=$1
BASE=$2
NUMBER=$3
SIGMA=$4
RAW=$5
INTERP=$6
if [ -z "$INTERP" ]; then
    INTERP=splineper
fi
BOUNDARY=$7
if [ -z "$BOUNDARY" ]; then
    BOUNDARY=hsym
fi

#maximal displacment of the image corners
L=$8
if [ -z "$L" ]; then
    L=3
fi

# type of transformation
# 2 --> translation
# 3 --> euclidean
# 6 --> affinity
# 8 --> homography
TYPE=$9
if [ -z "$TYPE" ]; then
    TYPE=8
fi

ZOOM=${10}
if [ -z "$ZOOM" ]; then
    ZOOM=1
fi

# create standard burst
create_burst $IN $BASE $NUMBER $INTERP $BOUNDARY $L $TYPE $ZOOM

# add noise
for i in `seq 1 $NUMBER`; do
    add_noise $SIGMA ${BASE}_$i.tiff ${BASE}_$i.tiff
done

# compute corresponding CFA images
if [ "$RAW" -eq "1" ]; then
    for i in `seq 1 $NUMBER`; do
        rgb2raw ${BASE}_$i.tiff ${BASE}_cfa_$i.tiff
    done
fi