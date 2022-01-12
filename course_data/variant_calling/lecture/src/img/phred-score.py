#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math


n = 1000
p = []
q = []

for i in range(n-1):
    p.append((i+1)*1.0/n)
    q.append(-10*math.log10(p[-1]))

fig = plt.figure(figsize=(3.5,2.5))
ax1 = fig.add_subplot(111)
ax1.plot(p,q, color='#d9534f',lw=2.5)
ax1.set_xlabel(r'$P_{err}$')
ax1.set_ylabel(r'$Q$')
ax1.set_xlim(0,1)
ax1.annotate(r'$Q = -10\,\log_{10}\, P_{err}$',xy=(0.4,0.7),xycoords='axes fraction')
plt.subplots_adjust(left=0.15,right=0.93,bottom=0.2,top=0.93)
plt.savefig('img/phred-score.pdf')


p = []
q = []
for i in range(n-1):
    p.append((i+1)*1.0/n)
    q.append(-10*math.log10(1-p[-1]))

fig = plt.figure(figsize=(3.5,2.5))
ax1 = fig.add_subplot(111)
ax1.plot(p,q, color='#d9534f',lw=2.5)
ax1.set_xlabel(r'$P(G|D)$')
ax1.set_ylabel(r'$Q$')
ax1.set_xlim(0,1)
ax1.annotate(r'$Q = -10\,\log_{10}\,[1 - P(G|D)]$',xy=(0.2,0.7),xycoords='axes fraction')
plt.subplots_adjust(left=0.15,right=0.93,bottom=0.2,top=0.93)
plt.savefig('img/phred-score-pdg.pdf')


plt.close()


