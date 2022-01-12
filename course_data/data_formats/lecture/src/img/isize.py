#!/usr/bin/env python
#
#   /nfs/users/nfs_p/pd3/sandbox/usr/x10-processing/NA12878
#   samtools stats /lustre/scratch113/projects/hiseqx_test/bwa_mem-hs37d5/13350/13350_1_mkdup.bam
#   cat rmme.bchk | grep ^IS | cut -f2,3 | gzip -c > isize.tab.gz
#   python isize.py
#

import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv,gzip

fname = 'isize.tab.gz'
dat = []
max  = 0
imax = 0

csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)
f = gzip.open(fname,'rb')
reader = csv.reader(f, 'tab')
for row in reader:
    x = float(row[0])
    y = float(row[1])
    dat.append([x,y])
    if max<y: max = y; imax = x

# shaded area
beg = imax*(1-0.25)
end = imax*(1+0.25)
verts = [(beg,0)]
for x in dat:
    if x[0]>=beg and x[0]<=end: verts.append(x)
verts.append((end,0))

fig, ax = plt.subplots(1, 1, figsize=(5,3.5))
ax.plot([x[0] for x in dat],[x[1] for x in dat],'-',color='k')
ax.plot([imax,imax],[0,max],'--',color='k')
ax.plot([verts[1][0],verts[1][0]],[0,verts[1][1]],'--',color='k')
ax.plot([verts[-2][0],verts[-2][0]],[0,verts[-2][1]],'--',color='k')
poly = Polygon(verts, facecolor='orange', edgecolor='orange')
ax.add_patch(poly)

ax.annotate('',(imax,1e6),(end,1e6),arrowprops={'arrowstyle':'<->'})
ax.annotate(r'1.25x',((imax+end)*0.5,1.2e6),ha='center')

ax.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
ax.set_yticks([])
ax.set_ylabel('Number of read pairs')
ax.set_xlabel('Insert Size')
ax.set_ylim(0,4.5e6)
ax.set_xlim(100)
ax.set_xticks([imax])
ax.set_xticklabels([r'x'])
ax.set_title('Read pairs within 25% of the main peak',fontsize=12)

plt.subplots_adjust(right=0.95,left=0.1,bottom=0.15)

plt.savefig('isize.pdf')
plt.close()

