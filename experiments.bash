#!/usr/bin/env bash

EDS_DIRECTORY="../vcf2eds/scripts/"
RESULTS_FILE="results.csv"
REPEATS=100

echo "filename,dataset,deg_prob,rec_prob,seg_size_avg,elm_len_avg,len,repeats,pattern_length,algorithm,total_time" > "$RESULTS_FILE"

# Run BNDM algorithm on all EDS data, using randomly generated pattern
EDS_FILE_PATTERN="synth-*.eds"
PATTERN_LENGTH=32
ALGORITHM="bndm"
find "$EDS_DIRECTORY" -type f -name "$EDS_FILE_PATTERN" -print0 | while read -d $'\0' FILE
do
    echo "Running BNDM-EDS: ./eds_search \"$FILE\" \"$REPEATS\" \"$PATTERN_LENGTH\" \"$ALGORITHM\"..."
    echo -n "\"$FILE\"," >> "$RESULTS_FILE"
    echo -n "$FILE" | sed -n 's/.*\/\([^\/]*\)_p=\([0-9.]*\)_r=\([0-9.]*\)_s=\([0-9.]*\)_e=\([0-9.]*\)_l=0*\([0-9]*\)\.eds/\1,\2,\3,\4,\5,\6,/p' >> "$RESULTS_FILE"
    echo -n "$REPEATS,$PATTERN_LENGTH,$ALGORITHM," >> "$RESULTS_FILE"
    ./eds_search "$FILE" "$REPEATS" "$PATTERN_LENGTH" "$ALGORITHM" | sed -n "s/Total time:\s*\([0-9.]*\) s/\1/p" >> "$RESULTS_FILE"
done

# Run BNDM-AA algorithm on DNA EDS data (skip protein EDS data), using static pattern
PROTEIN_FILE_PATTERN="synth-protein_*.eds"
PATTERN_LENGTH=10
PATTERN="QRGKGTCDQK" # Expands to 2 DNA patterns
ALGORITHM="bndm-aa"
find "$EDS_DIRECTORY" -type f -name "$EDS_FILE_PATTERN" ! -name "$PROTEIN_FILE_PATTERN" -print0 | while read -d $'\0' FILE
do
    echo "Running BNDM-EDS: ./eds_search \"$FILE\" \"$REPEATS\" \"$PATTERN_LENGTH\" \"$ALGORITHM\" \"$PATTERN\"..."
    echo -n "\"$FILE\"," >> "$RESULTS_FILE"
    echo -n "$FILE" | sed -n 's/.*\/\([^\/]*\)_p=\([0-9.]*\)_r=\([0-9.]*\)_s=\([0-9.]*\)_e=\([0-9.]*\)_l=0*\([0-9]*\)\.eds/\1,\2,\3,\4,\5,\6,/p' >> "$RESULTS_FILE"
    echo -n "$REPEATS,$PATTERN_LENGTH,$ALGORITHM," >> "$RESULTS_FILE"
    ./eds_search "$FILE" "$REPEATS" "$PATTERN_LENGTH" "$ALGORITHM" "$PATTERN" | sed -n "s/Total time:\s*\([0-9.]*\) s/\1/p" >> "$RESULTS_FILE"
done