#!/bin/bash

# Read input file line by line
while read -r INPUTNAME OUTPUTNAME IDENTIFIER; do
    FULL_INPUT_PATH="~/bruxhcc/XtoHH_MX-1000to4000/$INPUTNAME"
    FULL_OUTPUT_PATH="../data/XtoHH/MX-1000-4000/$OUTPUTNAME"
    TEST_OUTPUT_PATH="${FULL_OUTPUT_PATH/PAIReD/PAIReD_test}"
    
    echo "Processing: INPUT=$FULL_INPUT_PATH, OUTPUT=$FULL_OUTPUT_PATH, ID=$IDENTIFIER"
    
    python processFileToPAIReD.py "$FULL_INPUT_PATH" "$FULL_OUTPUT_PATH" --physicsprocess 25 -g Ellipse
    
    # Store the target directory and copy the test output to it
    TEST_DIR="../data/test/XtoHH/MX-1000-4000/"
    cp "$TEST_OUTPUT_PATH" "$TEST_DIR$OUTPUTNAME"
    
    echo "Finished processing $FULL_INPUT_PATH"
done < input_names/custom_inputs_XtoHH.txt
