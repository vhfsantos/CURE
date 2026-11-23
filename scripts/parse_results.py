#!/usr/bin/env python3
"""Convert uce_kit GFF matches into GeneRegion assignment lists."""

import re
import sys


def attributes(text):
    return dict(field.split('=', 1) for field in text.split(';') if '=' in field)


def transcript_genes(gff_path):
    result = {}
    if not gff_path:
        return result
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            info = attributes(fields[8])
            if fields[2].lower() in ('mrna', 'transcript', 'ncrna', 'lnc_rna',
                                     'rrna', 'trna', 'primary_transcript',
                                     'pseudogenic_transcript'):
                if info.get('ID') and info.get('Parent'):
                    result[info['ID']] = info['Parent'].split(',')[0]
    return result


def safe_id(value):
    return re.sub(r'[^A-Za-z0-9_.()\-]', '_', value)


def gene_id(info, transcripts):
    if info.get('gene') or info.get('gene_id'):
        return 'gene-' + safe_id(info.get('gene') or info['gene_id'])
    for parent in info.get('Parent', '').split(','):
        if parent in transcripts:
            return safe_id(transcripts[parent])
        if parent.startswith('gene-'):
            return safe_id(parent)
    dbxref = re.search(r'(?:^|,)GeneID:([^,;]+)', info.get('Dbxref', ''))
    if dbxref:
        return 'gene-' + safe_id(dbxref.group(1))
    raise ValueError('Cannot resolve exon gene from GFF attributes: ' + str(info))


def read_df(df_path, filt, output, region, gff_path=None):
    transcripts = transcript_genes(gff_path) if region == 'exon' else {}
    with open(df_path) as source, open(output, 'w') as destination:
        for line in source:
            row = line.rstrip('\n').split('\t')
            if not row or not row[0]:
                continue
            locus = row[0].split(filt)[0] + '.nexus' if filt else row[0] + '.nexus'
            if region == 'intergenic':
                destination.write(locus + '\n')
                continue
            info = attributes(row[4])
            if region == 'intron':
                gene = info.get('ID')
                if not gene:
                    raise ValueError('Intron without gene ID: ' + line)
                destination.write('\t'.join((locus, safe_id(gene))) + '\n')
                continue
            gene = gene_id(info, transcripts)
            exon = info.get('ID') or 'exon-' + '-'.join((info.get('Parent', 'unknown'),
                                                        row[2], row[3]))
            destination.write('\t'.join((locus, safe_id(exon), gene)) + '\n')


if __name__ == '__main__':
    read_df(*sys.argv[1:])
