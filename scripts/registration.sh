#!/bin/bash
# Registration of images
# For raw = 1 this is the two-step registration method for mosaicked images

# read input parameters
if [ "$#" -lt "4" ]; then
    echo "usage: $0 in_path homo_path number raw"
    echo "example: $0 im_%i.tiff homo_%i.hom 100 1"
    exit 1
fi

INPATH=$1
HOMOPATH=$2
NUMBER=$3
RAW=$4

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH

# Step 1: lowpass filtering
if [ "$RAW" -eq "1" ]; then
    echo "Applying the lowpass filter to the CFA images"

    INPATH2=registration_lowpass_%i.tiff
    for i in `seq 1 $NUMBER`; do
        INi=`printf $INPATH $i`
        OUTi=`printf $INPATH2 $i`
        echo "lowpass $INi $OUTi 0"
    done | parallel
else
    INPATH2=$INPATH
fi

# Step 2
echo "Starting standard registration"

ref_image=`printf $INPATH2 1`
for i in `seq 2 $NUMBER`; do
    INi=`printf $INPATH2 $i`
    REGi=`printf $HOMOPATH $i`
    to_echo="inverse_compositional_algorithm $ref_image $INi -f $REGi -o 1"
    if [ "$RAW" -eq "1" ]; then
        to_echo="$to_echo; rm $INi"
    fi
    echo "$to_echo"
done | parallel

if [ "$RAW" -eq "1" ]; then
    rm $ref_image
fi

