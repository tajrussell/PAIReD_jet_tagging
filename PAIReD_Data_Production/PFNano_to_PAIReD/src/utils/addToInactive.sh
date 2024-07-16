#!/bin/bash

# execute by shell command:
# $ bash addToTestset.sh $PATH $FILENO1 $FILENO2 $FILENO3 ...

# First argument is the path to the training set
# All following arguments are the index numbers of files to be added to the inactive list

echo "rename PAIReD files to inactive_PAIReD files"
# loop over all index numbers
for ind in "${@:2}"
do
    mv "$1/PAIReD_${ind}.root" "$1/inactive_PAIReD_${ind}.root"
done

echo "done"
