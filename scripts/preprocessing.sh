#!/bin/bash
# Preprocessing of the images
# The input filenames must be formatted as FILENAME%i.EXT where %i represents the number of the image
# The reference image must be the first one
# Set raw = 1 to handle RAW images
# For standard images raw_path is not used

# function for the ponomarenko noise estimation
ponomarenko_estimation() {
    base=$(echo ${file_pono_in%.*})
    split_raw $file_pono_in PNG16:$base png

    for j in 1 2 3 4; do
        ponomarenko ${base}_$j.png > ${file_pono_out}_$j.txt
        rm ${base}_$j.png
    done
}

# read input parameters
if [ "$#" -lt "6" ]; then
    echo "usage: $0 in_path raw_path im_path ind_ini ind_end raw crop_size"
    echo "example with RAW images: $0 PA%i.ORF raw_%i.tiff cfa_%i.tiff 1 100 1 512"
    echo "example with standard images: $0 PA%i.tiff /dev/null im_%i.tiff 1 100 0 512"
    exit 1
fi

INPATH=$1
RAWPATH=$2
IMPATH=$3
INI=$4
END=$5
RAW=$6
CROP=$7 # crop of the input images (set to 0 if no crop should be used)
if [ -z "$CROP" ]; then
    CROP=0
fi

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=${SCRIPTPATH%/*}/build/:$PATH

# number of input images
number=$(($END - $INI + 1))

# get the RAW images using DCRAW or just copy input
if [ "$RAW" -eq "1" ]; then
    ext=${INPATH##*.} # capture the file extension
    if [ "$ext" = "tiff" ]; then
        echo "Copying input RAW images"
        for i in `seq $INI $END`; do
            INi=`printf $INPATH $i`
            j=$(($i - $INI + 1))
            OUTj=`printf $RAWPATH $j`
            if [ ! -f $OUTj ]; then
                cp $INi $OUTj
            fi
        done
    else
        # transform from RAW to TIFF format
        echo "Transforming from RAW to TIFF format using DCRAW..."
        for i in `seq $INI $END`; do
            INi=`printf $INPATH $i`
            j=$(($i - $INI + 1))
            OUTj=`printf $RAWPATH $j`
            if [ ! -f $OUTj ]; then
                dcraw -v -4 -d -T -c $INi > $OUTj 2> /dev/null
            fi
        done
    fi
else
    RAWPATH=$IMPATH
    echo "Copying input images"
    for i in `seq $INI $END`; do
        INi=`printf $INPATH $i`
        j=$(($i - $INI + 1))
        OUTj=`printf $RAWPATH $j`
        if [ ! -f $OUTj ]; then
            cp $INi $OUTj
        fi
    done
fi

# crop the input images
# for standard images RAWPATH=IMPATH
if [ "$CROP" -gt "0" ]; then
    echo "Cropping the input images..."

        # compute the crop indices in order to keep the RGGB pattern
        RAWi=`printf $RAWPATH 1`
        w=$(identify -format "%w" $RAWi)
        h=$(identify -format "%h" $RAWi)
        iniw=$(($w/2 - 1 - $CROP/2))
        iniw=$(($iniw + $iniw%2))
        endw=$(($iniw + $CROP))
        inih=$(($h/2 - 1 - $CROP/2))
        inih=$(($inih + $inih%2))
        endh=$(($inih + $CROP))

        # perform the crop
        for i in `seq 1 $number`; do
            RAWi=`printf $RAWPATH $i`
            crop $iniw $inih $endw $endh $RAWi $RAWi
        done
fi

if [ "$RAW" -eq "1" ]; then
    # Transform of the images to correct the noise model (vst)
    echo "Applying the VST..."

    vst_type=0
    raw_ref=`printf $RAWPATH 1`
    file_pono_in=$raw_ref
    file_pono_out=pono
    ponomarenko_estimation

    # verify if there is enough points to evaluate the noise curve
    cmd="`cat ${file_pono_out}_1.txt`"
    l1=` echo "$cmd" | wc -l `
    if [ "$l1" -lt "3" ]; then
        echo "Not enough points to evaluate the noise curve"
        echo "Adding images..."
        
        # concatenate 4 images to have a larger image and a better noise curve estimation
        concatenated=concatenated.tiff
        raw2=`printf $RAWPATH 2`
        raw3=`printf $RAWPATH 3`
        raw4=`printf $RAWPATH 4`
        concatenate_images $raw_ref $raw2 $raw3 $raw4 $concatenated
        
        file_pono_in=$concatenated
        ponomarenko_estimation
        rm $concatenated
        
        cmd="`cat ${file_pono_out}_1.txt`"
        l1=` echo "$cmd" | wc -l `
        if [ "$l1" -lt "3" ]; then
            echo "Still not enough points to evaluate the noise curve..."
            rm ${file_pono_out}*
            exit 1
        fi
    fi
    RAWPATH2=${RAWPATH%%_*} # if RAWPATH=raw_%i.tiff then RAWPATH2=raw
    IMPATH2=${IMPATH%%_*} 
    ponomarenko_fit_raw_multiple $file_pono_out $RAWPATH2 tiff $IMPATH2 $vst_type 1 $number
    rm ${file_pono_out}*
fi

# Channel and mean equalization
echo "Channel and mean equalization"
    
    im_ref=`printf $IMPATH 1`
    if [ "$RAW" -eq "1" ]; then
        channel_equalization $im_ref $im_ref # multiplicative mean equalization of the color channels
    fi
    affine0255 $im_ref $im_ref 255 # set the max to 255

    eq_type=meanx
    for i in `seq 2 $number`; do
        IMi=`printf $IMPATH $i`
        equalization $im_ref $IMi $IMi $eq_type $RAW
    done

