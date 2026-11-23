#!/usr/bin/env python
# coding: utf-8
# Copy a GFF3 file to stdout and add introns between exons of the same transcript.

import sys
from collections import defaultdict


def attributes(field):
    return dict(part.split('=', 1) for part in field.split(';') if '=' in part)


def add_introns(infile, outfile, geneid_for_exon=False):
    transcript_gene = {}
    exons = defaultdict(list)
    current_gene = ''

    for line in infile:
        if not line or line.startswith('#'):
            outfile.write(line)
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            outfile.write(line)
            continue

        kind = fields[2].lower()
        attrs = attributes(fields[8])
        if kind == 'gene':
            current_gene = attrs.get('ID', '')
        elif kind in ('mrna', 'transcript', 'ncrna', 'lnc_rna', 'rrna',
                      'trna', 'primary_transcript', 'pseudogenic_transcript'):
            transcript_id = attrs.get('ID')
            gene_id = attrs.get('Parent', '').split(',')[0]
            if transcript_id and gene_id:
                transcript_gene[transcript_id] = gene_id
        elif kind == 'exon':
            if geneid_for_exon and current_gene:
                fields[8] = 'ID=' + current_gene
                line = '\t'.join(fields) + '\n'
            if fields[3].isdigit() and fields[4].isdigit():
                start, end = int(fields[3]), int(fields[4])
                if 0 < start <= end:
                    for parent in attrs.get('Parent', '').split(','):
                        if parent:
                            exons[(parent, fields[0], fields[6])].append((start, end))
        outfile.write(line)

    for (parent, scaffold, strand), intervals in exons.items():
        gene_id = transcript_gene.get(parent)
        if not gene_id and parent.startswith('gene-'):
            gene_id = parent
        if not gene_id:
            continue
        intervals.sort()
        last_end = intervals[0][1]
        for start, end in intervals[1:]:
            if start > last_end + 1:
                outfile.write('\t'.join((scaffold, 'Insert', 'intron',
                    str(last_end + 1), str(start - 1), '.', strand, '.',
                    'ID=' + gene_id)) + '\n')
            last_end = max(last_end, end)


def main(argv):
    gff_file = sys.stdin
    geneid_for_exon = False
    for arg in argv[1:]:
        if arg == '-h':
            sys.stderr.write('\n    add_introns_to_gff.py [<gff_file>] [-g]\n\n'
                             '    No <gff_file> reads from stdin.\n'
                             "    -g replaces exon descriptions with the gene's ID.\n\n")
            return
        elif arg == '-g':
            geneid_for_exon = True
        else:
            gff_file = open(arg, 'r')
    add_introns(gff_file, sys.stdout, geneid_for_exon)
    if gff_file is not sys.stdin:
        gff_file.close()


if __name__ == '__main__':
    main(sys.argv)
