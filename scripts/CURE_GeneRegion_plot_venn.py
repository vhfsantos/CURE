from venn import venn
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
ver = args[6]
locations = args[7]

def read_n_create_set(file):
    venn_list = list()
    with open(file, 'r') as handler:
        for loci in handler.readlines():
            venn_list.append(loci.strip())
    return set(venn_list)

def write_log(msg):
   sys.stdout.write("[ CURE v{} | ".format(ver) +
                 time.strftime("%Y-%m-%d %H:%M:%S") +
                 " ] " +
                 msg + '\n')

all_data = {'Assigned to intergenic regions': read_n_create_set(intergenic),
            'Assigned to exons': read_n_create_set(exons),
            'Assigned to introns': read_n_create_set(introns),
            'All UCEs in NEXUS dir': read_n_create_set(all_uces)}

# plot venn diagram
venn(all_data)
plt.savefig(os.path.join(output, 'CURE_stats.pdf'))

# Count only loci present in the input alignments. A mapped locus with an
# unsupported GFF feature (for example lncRNA) is neither unmapped nor
# confirmed intergenic.
input_uces = all_data['All UCEs in NEXUS dir']
exons = all_data['Assigned to exons'] & input_uces
introns = all_data['Assigned to introns'] & input_uces
intergenic = all_data['Assigned to intergenic regions'] & input_uces
with open(locations) as location_file:
    mapped = {line.split('\t', 1)[0] + '.nexus'
              for line in location_file if line.strip()} & input_uces
unclassified_mapped = mapped - exons - introns - intergenic
unassigned = input_uces - mapped

with open(os.path.join(output, 'CURE_stats.csv'),'w') as out_csv:
    out_csv.write('type, uce_count\n')
    out_csv.write('unassigned, {}\n'.format(len(unassigned)))
    out_csv.write('unclassified_mapped, {}\n'.format(len(unclassified_mapped)))
    out_csv.write('intron, {}\n'.format(len(introns - exons)))
    out_csv.write('exon, {}\n'.format(len(exons - introns)))
    out_csv.write('exon_and_intron, {}\n'.format(len(exons & introns)))
    out_csv.write('intergenic, {}\n'.format(len(intergenic)))

total_uces = len(input_uces)
write_log('------------------------------------------------------')
write_log('---------------------- Summary -----------------------')
write_log('------------------------------------------------------')
write_log('{} UCEs in NEXUS dir:'.format(str(total_uces)))
write_log('{} unassigned UCEs (not mapped to the reference)'.format(len(unassigned)))
write_log('{} mapped UCEs with no supported region type'.format(len(unclassified_mapped)))
write_log('{} assigned to introns'.format(len(introns - exons)))
write_log('{} assigned to exons'.format(len(exons - introns)))
write_log('{} assigned to exons and introns (accounted for exons)'.format(len(exons & introns)))
write_log('{} assigned to intergenic regions'.format(len(intergenic)))
write_log('------------------------------------------------------')
