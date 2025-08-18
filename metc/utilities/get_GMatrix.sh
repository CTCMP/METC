#!/bin/sh

echo "Please input filename: "
read input_file


output_file="gMatrix.dat"

awk '
/^\s*------------------------------------------------------------------------------\s*$/ {
    if (in_block) {
        print "";
        in_block = 0;
    } else {
        in_block = 1;
    }
    next;
}
in_block' "$input_file"  > "$output_file"

