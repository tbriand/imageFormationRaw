#!/bin/bash
# Evaluate the performance of the image formation algorithm
# A burst is built from an image and then the image is reconstructed
# Both registration and reconstruction precisions are evaluated

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH
PATH=${SCRIPTPATH%/*}/scripts/:$PATH

# read input parameters
if [ "$#" -lt "4" ]; then
    echo "usage: $0 in number sigma raw [interp boundary L type zoom]"
    echo "example: $0 im.tiff 100 5 0"
    echo "Specify the directory of temporary files using the variable DIR"
    echo "Use FREE=1 to free temporary files"
    exit 1
fi

IN=$1
NUMBER=$2
SIGMA=$3
RAW=$4
INTERP=$5
if [ -z "$INTERP" ]; then
    INTERP=p+s-spline11-spline3
fi
BOUNDARY=$6
if [ -z "$BOUNDARY" ]; then
    BOUNDARY=hsym
fi

# maximal displacment of the image corners
L=$7
if [ -z "$L" ]; then
    L=3
fi

# type of transformation
# 2 --> translation
# 3 --> euclidean
# 6 --> affinity
# 8 --> homography
TYPE=$8
if [ -z "$TYPE" ]; then
    TYPE=8
fi

ZOOM=$9
if [ -z "$ZOOM" ]; then
    ZOOM=1
fi

# make temporary directory for the experiments
if [ -z "$DIR" ]; then
    DIR=evaluate_performance_tmp
fi
mkdir $DIR

# remove the intermediate files to save memory (set to 1)
if [ -z "$FREE" ]; then
    FREE=0
fi

# create the burst
echo "Creating synthetic burst"
base=$DIR/burst
create_burst.sh $IN $base $NUMBER $SIGMA $RAW $INTERP $BOUNDARY $L $TYPE $ZOOM

# registration
echo "Starting registration step"
homo_path=$DIR/estimated_%i.hom 
if [ "$RAW" -eq "1" ]; then
    im_path=${base}_cfa_%i.tiff
else
    im_path=${base}_%i.tiff
fi
registration.sh $im_path $homo_path $NUMBER $RAW

# image fusion
echo "Starting fusion step" 
out=$DIR/reconstructed.tiff
fusion_irregularly_sampled_data.sh $im_path $homo_path $out $NUMBER $RAW $ZOOM

# evaluation of the registration precision
error=$DIR/registration_error.txt
rm -f $error
w=`identify -format %w $IN`
h=`identify -format %h $IN`
field=$DIR/field.tiff
for i in `seq 2 $NUMBER`; do
    estimated=`printf $homo_path $i`
    truth=${base}_$i.hom
    compare_homography $w $h "`cat $estimated`" "`cat $truth`" $field 1
    compute mean 0 $field >> $error
done
rm $field
echo "Precision of the registration (EPE):"
mean_and_std $error 0
        
# evaluation of the reconstruction precision
cropv=20
ref=$DIR/ref.tiff
crop $cropv $cropv -$cropv -$cropv $IN $ref
diff=$DIR/diff.tiff
difference_images $out $ref $diff
echo "Precision of the reconstruction (RMSE):"
compute rmse 0 $diff

# if asked, free the intermediate files
if [ "$FREE" -eq "1" ]; then
    rm -r $DIR
fi