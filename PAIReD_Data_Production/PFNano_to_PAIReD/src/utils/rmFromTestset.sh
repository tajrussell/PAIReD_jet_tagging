#!/bin/bash

# execute by shell command:
# $ bash rmFromTestset.sh $PATH $FILENO1 $FILENO2 $FILENO3 ...

# First argument is the path to the training set
# All following arguments are the index numbers of files to be removed to the test set

echo "rename PAIReD files to test_PAIReD files"
# loop over all index numbers
for ind in "${@:2}"
do
    mv "$1/test_PAIReD_${ind}.root" "$1/PAIReD_${ind}.root"
done

echo "done"
