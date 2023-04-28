#!/usr/bin/env python
# coding: utf-8
# read from stdin and write to stdout inserting introns between exon lines of a gff3 format file
# change to read from stdin if no args or assume first arg is the gff file on which to add introns
# use -g option to replace the exon line descriptions with the appropriate geneid

import sys

def add_introns(infile, outfile, geneid_for_exon = False):
    def to_int(digits_str):  # do not throw a ValueError exception if input happens to not be a string of just digits
        return int(digits_str) if digits_str.isdigit() else 0

    cur_flds = []; geneid = ""
    for gff_ln in infile:
        if len(gff_ln) < 3 or gff_ln[0] == "#":  # output empty or comment lines
            outfile.write(gff_ln)
            continue
        
        lst_flds = list(cur_flds)  # copy the last non-comment line fields
        cur_flds = gff_ln.split("\t")
        
        too_short = len(cur_flds) < 5 or len(lst_flds) < 5
        if too_short:
            outfile.write(gff_ln)
            continue
        
        cur_type = cur_flds[2].lower()
        if cur_type == "gene":  # remember the geneid to associate it with the intron
            geneid = cur_flds[8].split(";")[0]
        
        if cur_type != "exon" or lst_flds[2].lower() != "exon":
            if geneid_for_exon and cur_type == "exon":
                cur_flds[8] = geneid
                gff_ln = "\t".join(cur_flds) + "\n"
                
            outfile.write(gff_ln)
            continue
        
        # last line and current line both are exon lines, this is where to insert the intron line        
        cur_beg = to_int(cur_flds[3]); cur_end = to_int(cur_flds[4])
        lst_beg = to_int(lst_flds[3]); lst_end = to_int(lst_flds[4])
        
        if cur_beg==0 or cur_end==0 or lst_beg==0 or lst_end==0:
            # we expected positive integers but didn't get them
            outfile.write(gff_ln)
            continue
        
        if lst_end < cur_beg:  # positive strand
            int_beg = lst_end + 1
            int_end = cur_beg - 1
        else:  # negative strand
            int_beg = cur_end + 1
            int_end = lst_beg - 1
        
        intron = "\t".join([ lst_flds[0], "Insert","intron", str(int_beg),str(int_end), lst_flds[5],lst_flds[6],lst_flds[7], geneid ]) + "\n"
        outfile.write(intron)  # output intron line
        
        if geneid_for_exon and cur_type == "exon":
             cur_flds[8] = geneid
             gff_ln = "\t".join(cur_flds) + "\n"

        outfile.write(gff_ln)  # output current exon line
        
def main(argv):
    gff_file = sys.stdin
    geneid_for_exon = False
    for ix in range(1, len(argv)):
        arg = argv[ix]
        if arg == "-h":
            usage = "\n    add_introns_to_gff.py [<gff_file>] [-g]\n\n    No <gff_file> reads from stdin\n    -g will replace exon description with gene's ID\n\n"
            sys.stderr.write(usage)
            sys.exit(0)
        elif arg == "-g":
            geneid_for_exon = True
        else:
            gff_file = arg

    if gff_file == sys.stdin:
        sys.stderr.write("reading gff from stdin. for usage, use -h option.\n")
    else:
        gff_file = open(gff_file, 'r')

    add_introns(gff_file, sys.stdout, geneid_for_exon)

if __name__ == '__main__':
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) # don't complain if pipe closes output (head or less commands will cause this)
    
    try:
        main(sys.argv)
    except (KeyboardInterrupt, SystemExit): # don't do stack trace for Ctrl+C
        sys.exit(0)
