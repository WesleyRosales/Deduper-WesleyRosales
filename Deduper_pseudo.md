## Define the Problem
PCR duplicates are sequences that originate from the same DNA molecule. This introduces introduces false sequences in the dataset and must be removed prior to analysis.

## Examples
```
14 test input rows (2 for each category) that test:                             Expected Action:
Soft clipping adjusted to match another non-soft clipped sequence               Keep 1 (the non-soft clipped sequence)
Soft clipping that adjust to unique rows                                        Keep 2
Same position and different UMI                                                 Keep 2
Same UMI and different position                                                 Keep 2
Same UMI and same position                                                      Keep 1 (first one)
Different UMI and different position                                            Keep 2
Error in UMI and different position                                             Keep 0
                                                                                Total: 10 rows in output
```

## Pseudocode
PCR duplicates can be identified by having the same chromosome (col 3), position (col 4), strand (col 2), and Unique Molecular Index (UMI). This program should check for matches in these 4 categories as well as consider soft clipping in the algorithm.

1. Use Samtools to sort by position
2. Argparse for input .sam file, UMI file, and output filename.
3. Write all header rows to output file.
4. Save alignments with the same position until a row with a new position is found in "alignments" object. Save the position and UMI of contents of list (chromosome, position, strand, UMI) to another list that remains until end of deduplication named position list.
5. Check for any soft clipping. If found, adjust position and check if position + UMI is in position list -> if not in position list, write to output file and save position+UMI in position list. Remove soft clipping from list.
7. Divide "alignments" by strand and run steps 6 and 7 in parallel.
6. Filter out alignments with same UMIs AND with UMIs NOT IN the UMI input file.
7. Write remaining alignments to input and remove from list of lists
8. Repeat process for next position.
9. Once end of file is reached, close input and output files.

## High-level Functions
Name: findUMI <br>
Function: Isolate UMI from header column in sam file <br>
Test Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC, NS500451:154:HWKTMBGXX:1:11101:18996:1145:TTCGCCTA <br>
Test Output: CTGTTCAC, TTCGCCTA