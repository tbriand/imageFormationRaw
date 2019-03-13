#!/bin/bash
# Fusion of irregularly sampled data by accumulation

# read input parameters
if [ "$#" -lt "4" ]; then
    echo "usage: $0 inpath homo_path output_image number zoom"
    echo "example: $0 cfa_%i.tiff homo_%i.tiff out.tiff 101 2"
    exit 1
fi

INPATH=$1
HOMOPATH=$2
OUT=$3
NUMBER=$4
ZOOM=$5
if [ -z "$ZOOM" ]; then
    ZOOM=1
fi

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH

# set parameters
ind_ini=2 # change to 1 to use the reference image during the fusion
sigma=0.70710678118 # scale of the gaussian weight (sqrt(2)/2)
order=0
if [ "$ZOOM" -gt "1" ]; then
    order=2
fi

# irregularly sampled data fitting
echo "Starting data accumulation and blurry image computation"

    INPATH2=${INPATH%%_*}
    HOMOPATH2=${HOMOPATH%%_*}
    combi_ckr $INPATH2 $HOMOPATH2 $ind_ini $NUMBER $ZOOM $OUT order $sigma 0

echo "Done"

# sharpening step
echo "Sharpening step"

    cropv=10 # the image is cropped before and after the sharpening step to remove the boundary artefacts
    crop $cropv $cropv -$cropv -$cropv $OUT $OUT
    asymptotic_nc_filter $OUT $OUT $order $sigma 1
    crop $cropv $cropv -$cropv -$cropv $OUT $OUT

echo "Done"
