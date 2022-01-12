import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

sample_ids = False

dat = []
#with open('gtcheck.fail/gtcheck.tab', 'rb') as f:
#with open('gtcheck1/gtcheck.tab', 'rb') as f:
#with open('gtcheck2/gtcheck.tab', 'rb') as f:
with open('gtcheck3/gtcheck.tab', 'rb') as f:
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0][0]=='#': continue
        if row[0]!='CN': continue
        tgt = 0
        if row[4]=='HG00096': tgt = 1
        dat.append([float(row[1]), float(row[2]), float(row[3]), tgt, row[4]])

dat = sorted(dat)

iq = -1; dp = 0
for i in range(len(dat)):
    if iq==-1 and dat[i][3]==1: iq = i
    dp += dat[i][2]
dp /= len(dat)

fig,ax1 = plt.subplots(figsize=(6,4))
ax2 = ax1.twinx()
plots  = ax1.plot([x[0] for x in dat],'o-', ms=5, color='g', mec='g', label='Discordance (total)')
plots += ax1.plot([x[1] for x in dat], '^', ms=5, color='r', mec='r', label='Discordance (avg per site)')
plots += ax2.plot([x[2] for x in dat],'v', ms=5, color='k', label='Number of sites')
if iq!=-1:
   plots += ax1.plot([iq],[dat[iq][0]],'o',color='orange', ms=10, label='Tested sample')
for tl in ax1.get_yticklabels(): tl.set_color('g')
for tl in ax2.get_yticklabels(): tl.set_color('k'); tl.set_fontsize(9)
min_dp = min([x[2] for x in dat])
max_dp = max([x[2] for x in dat])
ax2.set_ylim(min_dp-1,max_dp+1)
ax1.set_xlim(-0.05*len(dat),1.05*(len(dat)-1))
ax1.set_xlabel('Sample ID')
plt.subplots_adjust(left=0.1,right=0.9,bottom=0.12,top=0.9)
if sample_ids:
   ax1.set_xticks(range(len(dat)))
   ax1.set_xticklabels([x[4] for x in dat],**{'rotation':45, 'ha':'right', 'fontsize':8})
   plt.subplots_adjust(bottom=0.2)
ax1.set_ylabel('Discordance',color='g')
ax2.set_ylabel('Number of sites',color='k')
ax2.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')
labels = [l.get_label() for l in plots]
#plt.legend(plots,labels,numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
plt.savefig('gtcheck1/gtcheck.png')
plt.savefig('gtcheck1/gtcheck.pdf')
plt.close()

