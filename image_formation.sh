#!/bin/bash
# Main script for the image formation (from RAW images or not)
# The input filenames must be formatted as FILENAME%i.EXT where %i represents the number of the image
# The reference image must be the first one
# Set raw = 1 to handle RAW images

# read input parameters
if [ "$#" -lt "5" ]; then
    echo "usage: $0 in_path ind_ini ind_end output_image raw zoom crop_size free"
    echo "example: $0 PA%i.ORF 1 100 out.tiff 1 1 512 0"
    exit 1
fi

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=$SCRIPTPATH/scripts/:$PATH
PATH=$SCRIPTPATH/build/:$PATH

# required parameters
INPATH=$1
INI=$2
END=$3
OUT=$4
RAW=$5

# optional parameters
ZOOM=$6
if [ -z "$ZOOM" ]; then
    ZOOM=1
fi
CROP=$7 # crop of the input images (set to 0 if no crop should be used)
if [ -z "$CROP" ]; then
    CROP=0
fi
FREE=$8 # remove the intermediate files to save memory (set to 1)
if [ -z "$FREE" ]; then
    FREE=0
fi

# number of input images
number=$(($END - $INI + 1)) 

# intermediate files are saved in the tmp directory
dir=tmp
mkdir $dir 
im_path=$dir/im_%i.tiff
homo_path=$dir/homography_%i.hom
raw_path=$dir/raw_%i.tiff

if [ "$RAW" -eq "1" ]; then
    echo "Starting image formation from RAW images"
else
    echo "Starting image formation from standard images"
fi

# preprocessing
echo "Starting preprocessing step"
preprocessing.sh $INPATH $raw_path $im_path $INI $END $RAW $CROP

# registration
echo "Starting registration step"
registration.sh $im_path $homo_path $number $RAW

# image fusion
echo "Starting fusion step"
fusion_irregularly_sampled_data.sh $im_path $homo_path $OUT $number $RAW $ZOOM

# if asked, free the intermediate files
if [ "$FREE" -eq "1" ]; then
    rm -r $dir
fi
