#!/bin/bash
# Two-step registration method for mosaicked images

# read input parameters
if [ "$#" -lt "3" ]; then
    echo "usage: $0 in_path homo_path number"
    echo "example: $0 cfa_%i.tiff homo_%i.hom 100"
    exit 1
fi

INPATH=$1
HOMOPATH=$2
NUMBER=$3

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH

# Step 1: lowpass filtering
echo "Applying the lowpass filter to the CFA images"

    low_path=registration_lowpass_%i.tiff
    for i in `seq 1 $NUMBER`; do
        INi=`printf $INPATH $i`
        OUTi=`printf $low_path $i`
        lowpass $INi $OUTi 0
    done

echo "Done"

echo "Starting grayscale registration"

    ref_image=`printf $low_path 1`
    for i in `seq 2 $NUMBER`; do
        INi=`printf $low_path $i`
        REGi=`printf $HOMOPATH $i`
        inverse_compositional_algorithm $ref_image $INi -f $REGi -o 1
        rm $INi
    done
    rm $ref_image

echo "Done"

