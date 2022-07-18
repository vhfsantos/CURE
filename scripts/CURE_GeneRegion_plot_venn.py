from venn import venn, generate_petal_labels, draw_venn, generate_colors
import matplotlib.pyplot as plt
import sys
import os
import time

# parse results. all below are files
args = sys.argv
intergenic = args[1]
exons = args[2]
introns = args[3]
all_uces = args[4]
output = args[5]

def read_n_create_set(file):
    venn_list = list()
    with open(file, 'r') as handler:
        for loci in handler.readlines():
            venn_list.append(loci.strip())
    return set(venn_list)

def write_log(msg):
   out_csv.write("[ CURE v{} | " + 
                 time.strftime("%Y-%m-%d %H:%M:%S") + 
                 " ]".format(version) +
                 msg)

all_data = {'Assigned to intergenic regions': read_n_create_set(intergenic),
            'Assigned to exons': read_n_create_set(exons),
            'Assigned to introns': read_n_create_set(introns),
            'All UCEs': read_n_create_set(all_uces)}

# plot venn diagram
venn(all_data)
plt.savefig(os.path.join(output, 'CURE_stats.pdf')) 

# saving stats in csv file
lbs = generate_petal_labels(all_data.values(), fmt="{size}")

with open(os.path.join(output, 'CURE_stats.csv'),'w') as out_csv:
    write_log('type, uce_count')
    write_log('unassigned, {}'.format(lbs['0001']))
    write_log('intron, {}'.format(lbs['0011']))
    write_log('exon, {}'.format(lbs['0101']))
    write_log('exon_and_intron, {}'.format(lbs['0111']))
    write_log('intergenic, {}'.format(lbs['1001']))

sys.stdout.write("""
Unassigned UCEs: {}\n
Assigned to introns: {}\n
Assigned to exons: {}\n
Assigned to exons and introns: {}\n
Assigned to intergenic regions: {}\n
""".format(lbs['0001'], lbs['0011'], lbs['0101'], lbs['0111'], lbs['1001']))