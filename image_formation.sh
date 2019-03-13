#!/bin/bash
# Main script for the image formation from RAW images
# The input filenames must be formatted as FILENAME%i.EXT where %i represents the number of the image
# The reference image must be the first one

# read input parameters
if [ "$#" -lt "4" ]; then
    echo "usage: $0 in_path ind_ini ind_end output_image zoom crop_size free"
    echo "example: $0 PA%i.ORF 1 100 out.tiff 1 512 0"
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

# optional parameters
ZOOM=$5
if [ -z "$ZOOM" ]; then
    ZOOM=1
fi
CROP=$6 # crop of the input raw images (set to 0 if no crop should be used)
if [ -z "$CROP" ]; then
    CROP=0
fi
FREE=$7 # remove the intermediate files to save memory (set to 1)
if [ -z "$FREE" ]; then
    FREE=0
fi

# number of input images
number=$(($END - $INI + 1)) 

# intermediate files are saved in the tmp directory
dir=tmp
mkdir $dir 
raw_path=$dir/raw_%i.tiff
cfa_path=$dir/cfa_%i.tiff
homo_path=$dir/homography_%i.hom

# preprocessing
echo "Starting preprocessing step"
preprocessing_raw_images.sh $INPATH $raw_path $cfa_path $INI $END $CROP

# registration
echo "Starting registration step"
two_step_registration.sh $cfa_path $homo_path $number

# image fusion
echo "Starting fusion step"
fusion_irregularly_sampled_data.sh $cfa_path $homo_path $OUT $number $ZOOM

# if asked, free the intermediate files
if [ "$FREE" -eq "1" ]; then
    rm -r $dir
fi
