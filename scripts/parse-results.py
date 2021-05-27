#!/usr/bin/env python
import sys

def read_df(df_path, filt, output, region):
    # open output and input files
    with open(df_path,'r') as df:
        with open(output, 'w') as outdf:
            for line in df.readlines():
                # input is tab delimited:
                row = line.strip().split('\t')
                # all types will get the UCE name filtered
                UCE_name = row[0]
                UCE_name_filtered = UCE_name.split(filt)[0]+'.nexus'
                # now, each type will have diff. treatments
                if region == "exon":
                    full_ID = row[4]
                    # full ID example for exons:
                    # BASE:GB51190;gbkey=mRNA;gene=LOC411181;product=sodium...
                    # GeneID is between 'gene=' and ';'
                    first_split = full_ID.split('gene=')[1]
                    second_split = first_split.split(';')[0]
                    gene_ID = 'gene-'+second_split
                    # same for exon ID:
                    # ID=exon-XM_016912896.2-10;Parent=rna-XM_0169128...
                    # ExonID is between 'ID=' and ';'
                    first_split = full_ID.split('ID=')[1]
                    exon_ID = first_split.split(';')[0]
                    # all done, let's write it
                    line2write = '\t'.join([UCE_name_filtered, exon_ID, gene_ID])+'\n'
                    outdf.write(line2write)
                elif region == "intron":
                    full_ID = row[4]
                    # full ID example for introns:
                    # ID=gene-LOC410882
                    # Just need to remove 'ID='
                    gene_ID = full_ID.split('ID=')[1]
                    # all done, let's write it
                    line2write = '\t'.join([UCE_name_filtered, gene_ID])+'\n'
                    outdf.write(line2write)
                elif region == 'intergenic':
                    # nothing to do, just write UCE name filtered
                    outdf.write(UCE_name_filtered+'\n')

def main():
    # read arguments
    args = sys.argv
    df_path = args[1]
    filt = args[2]
    output = args[3]
    region = args[4]
    print('df: {}'. format(df_path))
    print('filt: {}'. format(filt))
    print('output: {}'. format(output))
    print('region: {}'. format(region))
    # call function
    read_df(df_path, filt, output, region)
    
if __name__ == '__main__':
    main()