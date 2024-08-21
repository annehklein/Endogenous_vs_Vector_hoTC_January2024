#!/usr/bin/env python3

# this scirpt filters reads by a barcode that is present at the start of the read
# it also checks for the reverse complement at the end of the read
# it outputs reads with the barcode to the output file
# it outputs counts of total reads and kept reads into an optional log file

import argparse
import sys
import pdb

def main():
    
    # command line arguments
    parser = argparse.ArgumentParser(description = "Filter a fastq by barcode")
    parser.add_argument('-i', help='Input fastq', required=True)
    parser.add_argument('-o', help='Output fastq', required=True)
    parser.add_argument('-seq', help='Barcode sequence', required=True)
    parser.add_argument('-log', help='Log file')
    
    args = parser.parse_args()
    
    # open input and output fastq
    with open(args.i, 'r') as infile, open(args.o, 'w') as outfile:
        # to keep count of total and kept reads
        count = 0
        total = 0
        # buffer to store lines from one read
        buffer = []
        
        # iterate over input file
        for i, line in enumerate(infile):
            
            # if we've reached a new read
            if i % 4 == 0:
                # check if we want to keep it
                count += check_buffer(args.seq, buffer, outfile)
                # reset buffer
                buffer = []
                total += 1
            
            # add line to buffer
            buffer.append(line)
            
    # output info about how many reads were seen and kept
    if args.log:
        with open(args.log, 'w') as log:
            log.write(f"found {args.seq} in {count} of {total} reads ({count/total*100}%)\n")
            
            
def check_buffer(seq, buffer, outfile):
    """
    Check if the read in the buffer contains the barcode (seq)
    If it does, write it to outfile
    """
    
    if len(buffer) == 0:
        return 0
    
    # check start of read for barcode
    if buffer[1][:len(seq)].upper() == seq.upper():
        # trim barcode from sequence and qualities
        buffer[1] = buffer[1][len(seq):]
        buffer[3] = buffer[3][len(seq):]
        # write lines
        outfile.writelines(buffer)
        return 1
        
    # check end of read for reverse complement
    elif buffer[1].strip("\n")[-len(seq):].upper() == reverse_complement(seq):
        # trim barcode from end
        buffer[1] = buffer[1].strip("\n")[:-len(seq)] + "\n"
        buffer[3] = buffer[3].strip("\n")[:-len(seq)] + "\n"
        # write lines
        outfile.writelines(buffer)
        return 1
        
    return 0
    
    
def reverse_complement(seq):
    """
    Convert to upper case and reverse complement
    """
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars, replace_chars)
    
    return seq.upper().translate(tab)[::-1]

if __name__ == "__main__":
    main()

