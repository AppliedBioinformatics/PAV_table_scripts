# python makePAV.py input_bed_reference input_formatted_mosdepth_file

import sys
from collections import defaultdict
gene_dict = defaultdict(int) # key: gene name, value: total length of gene based on exons
with open(sys.argv[1]) as fh:
   for line in fh:
       ll = line.split()
       if len(ll) < 3: continue
       if ll[7] != 'CDS': continue
       start, end = map(int, ll[1:3])
       length = end - start
       name = ll[-1].split(';')[-2]
       assert 'Parent' in name
       name = name.replace('Parent=',)
       gene_dict[name] += length
f = sys.argv[2]
overlap_dict = defaultdict(int)
with open(f) as fh:
   for line in fh:
       ll = line.split()
       cov = int(ll[3])
       if cov < 2: continue
       start_overlap, end_overlap = map(int, ll[1:3])
       length_overlap = end_overlap - start_overlap
       name = ll[-2].split(';')[-2]
       assert 'Parent' in name
       name = name.replace('Parent=',)
       overlap_dict[name] += length_overlap
genes = sorted(gene_dict)
header = ['Individual'] + genes
print('\t'.join(header))
ll = [f]
for g in genes:
   overlap_percentage = overlap_dict[g]/float(gene_dict[g])*100
   if overlap_percentage <= 5:
       # LOST
       ll.append('0')
   else:
       ll.append('1')
print('\t'.join(ll))
