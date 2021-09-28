#!/usr/bin/env bash

EDS_DIRECTORY="../vcf2eds/scripts/"
EDS_FILE_PATTERN="synth-*.eds"
PROTEIN_FILE_PATTERN="synth-protein_*.eds"
NODEG_FILE_PATTERN="synth-nodeg*.eds"
SEG_SIZE_FILE_PATTERN="synth-seg-size-avg*.eds"
RESULTS_FILE="results_v2.csv"
REPEATS=10

echo "filename,dataset,deg_prob,rec_prob,seg_size_avg,elm_len_avg,len,repeats,pattern,pattern_length,algorithm,total_time" > "$RESULTS_FILE"

# Run BNDM algorithm on all EDS data (DNA and protein), using randomly generated pattern
PATTERN_LENGTH=30
ALGORITHM="bndm"
find "$EDS_DIRECTORY" -type f -name "$EDS_FILE_PATTERN" ! -name "$NODEG_FILE_PATTERN" -print0 | while read -d $'\0' FILE
do
    echo "Running BNDM-EDS: ./eds_search \"$FILE\" \"$REPEATS\" \"$PATTERN_LENGTH\" \"$ALGORITHM\"..."
    echo -n "\"$FILE\"," >> "$RESULTS_FILE"
    echo -n "$FILE" | sed -n 's/.*\/\([^\/]*\)_p=\([0-9.]*\)_r=\([0-9.]*\)_s=\([0-9.]*\)_e=\([0-9.]*\)_l=0*\([0-9]*\)\.eds/\1,\2,\3,\4,\5,\6,/p' >> "$RESULTS_FILE"
    echo -n "$REPEATS,,$PATTERN_LENGTH,$ALGORITHM," >> "$RESULTS_FILE"
    ./eds_search "$FILE" "$REPEATS" "$PATTERN_LENGTH" "$ALGORITHM" | sed -n "s/Total time:\s*\([0-9.]*\) s/\1/p" >> "$RESULTS_FILE"
done

# Run BNDM-AA algorithm on DNA EDS data (skip protein EDS data)
ALGORITHM="bndm-aa"
find "$EDS_DIRECTORY" -type f -name "$EDS_FILE_PATTERN" ! -name "$PROTEIN_FILE_PATTERN" ! -name "$NODEG_FILE_PATTERN" -print0 | while read -d $'\0' FILE
do
  while read PATTERN; do
    PATTERN_LENGTH=${#PATTERN}
    echo "Running BNDM-EDS: ./eds_search \"$FILE\" \"$REPEATS\" \"$PATTERN_LENGTH\" \"$ALGORITHM\" \"$PATTERN\"..."
    echo -n "\"$FILE\"," >> "$RESULTS_FILE"
    echo -n "$FILE" | sed -n 's/.*\/\([^\/]*\)_p=\([0-9.]*\)_r=\([0-9.]*\)_s=\([0-9.]*\)_e=\([0-9.]*\)_l=0*\([0-9]*\)\.eds/\1,\2,\3,\4,\5,\6,/p' >> "$RESULTS_FILE"
    echo -n "$REPEATS,$PATTERN,$PATTERN_LENGTH,$ALGORITHM," >> "$RESULTS_FILE"
    ./eds_search "$FILE" "$REPEATS" "$PATTERN_LENGTH" "$ALGORITHM" "$PATTERN" | sed -n "s/Total time:\s*\([0-9.]*\) s/\1/p" >> "$RESULTS_FILE"
  done < patterns_aa.list
done

# Run BNDM-AA and BNDM-EDS on non-degenerate data - to compare BNDM performance, using fixed DNA pattern (the same
# pattern is used as both patterns in BNDM-AA)
PATTERN="TTTAATCATGAGCTGTTGCATACGGCGC" # Single pattern to be used
ALGORITHMS="bndm bndm-aa"
find "$EDS_DIRECTORY" -type f -name "$NODEG_FILE_PATTERN" -print0 | while read -d $'\0' FILE
do
  for ALGORITHM in $ALGORITHMS; do
    PATTERN_LENGTH=${#PATTERN}
    echo "Running BNDM-EDS: ./eds_search \"$FILE\" \"$REPEATS\" \"$PATTERN_LENGTH\" \"$ALGORITHM\" \"$PATTERN\" \"$PATTERN\"..."
    echo -n "\"$FILE\"," >> "$RESULTS_FILE"
    echo -n "$FILE" | sed -n 's/.*\/\([^\/]*\)_p=\([0-9.]*\)_r=\([0-9.]*\)_s=\([0-9.]*\)_e=\([0-9.]*\)_l=0*\([0-9]*\)\.eds/\1,\2,\3,\4,\5,\6,/p' >> "$RESULTS_FILE"
    echo -n "$REPEATS,$PATTERN,$PATTERN_LENGTH,$ALGORITHM," >> "$RESULTS_FILE"
    ./eds_search "$FILE" "$REPEATS" "$PATTERN_LENGTH" "$ALGORITHM" "$PATTERN" "$PATTERN" | sed -n "s/Total time:\s*\([0-9.]*\) s/\1/p" >> "$RESULTS_FILE"
  done
done