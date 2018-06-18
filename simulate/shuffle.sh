#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
output_file=""
left_file=""
right_file=""
verbose=0

while getopts "h?l:r:" opt; do
    case "$opt" in
	h|\?)
	    echo "help"
	    exit 0
	    ;;
	l)  left_file=$OPTARG
	    ;;
	r)  right_file=$OPTARG
	    ;;
    esac
done

shift $((OPTIND-1))

base_left=`basename $left_file`
base_right=`basename $right_file`
dirname_left=`dirname $left_file`
dirname_right=`dirname $right_file`

if [ $dirname_left == $dirname_right ]
then
    fl=${base_left%*.fasta}_shuffled.fa
    fr=${base_right%*.fasta}_shuffled.fa
    out_left="$dirname_left/${fl}"
    out_right="$dirname_right/${fr}"
    cmd="/usr/bin/time paste ${left_file} ${right_file} | paste - - | shuf | awk '{print \$1 > \"$out_left\" ; print \$3 > \"$out_left\"; print \$2 > \"$out_right\"; print \$4 > \"$out_right\"}'"
    echo $cmd
    eval $cmd
else
    echo "The left and right reads should be in the same directory"
    exit 1
fi
