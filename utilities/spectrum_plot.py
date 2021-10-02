#!/usr/bin/env python

from pylab import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

matplotlib.rcParams['lines.linewidth']=0.5
# matplotlib.rcParams['figure.dpi'] = 500
# matplotlib.rcParams['savefig.dpi'] = 500
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

sdata=np.loadtxt(sys.argv[1])

xval = np.zeros(sdata.shape[0], dtype=float)
yval = np.zeros((sdata.shape[1]-1,sdata.shape[0]), dtype=float)

for ipoint in range(sdata.shape[0]):
  xval[ipoint] = sdata[ipoint,0]
  for iorb in range(sdata.shape[1]-1):
    yval[iorb][ipoint] = sdata[ipoint,iorb+1]

plt.title("spectrum function")
plt.xlabel(r'$\omega$(eV)')
plt.ylabel("")
plt.xlim(-10.0,10.0)
plt.xticks(np.arange(-10.0,10.1, step=2.0))
spectrum = plt.subplot(1,1,1)
for iorb in range(sdata.shape[1]-1):
  spectrum.plot(xval,yval[iorb])

plt.savefig('spectrum_aims.png')
plt.show()

