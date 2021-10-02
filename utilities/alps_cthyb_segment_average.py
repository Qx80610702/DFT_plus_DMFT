#!/usr/bin/env python
'''
If the system is non-magnetic or paramagnetic, calculte average of 
spin-up and spin-down of self-energy, Matsubara function and imaginary
time Green function 
'''
import numpy as np

#====================
#     slef-energy
#====================
sdata=np.loadtxt("Sw.dat")
nomega=sdata.shape[0]
for iomega in range(nomega):
    m_tot=int(sdata.shape[1]-1)/4)
    for m in range(m_tot):
      sdata[iomega, 1+4*m] = (sdata[iomega, 1+4*m] + sdata[iomega, 3+4*m])/2.0
      sdata[iomega, 3+4*m] = sdata[iomega, 1+4*m]
      sdata[iomega, 2+4*m] = (sdata[iomega, 2+4*m] + sdata[iomega, 4+4*m])/2.0
      sdata[iomega, 4+4*m] = sdata[iomega, 2+4*m]

fi=open("Sw.out",'w')
for iomega in range(nomega):
    m_tot=int(sdata.shape[1]-1)/4)
    print("{:22.15f}".format(sdata[iomega, 0]), file=fi, end='')
    for m in range(m_tot):
        print("{:20.15f}".format(sdata[iomega, 1+4*m]), file=fi, end='')
        print("{:20.15f}".format(sdata[iomega, 2+4*m]), file=fi, end='')
        print("{:20.15f}".format(sdata[iomega, 3+4*m]), file=fi, end='')
        print("{:20.15f}".format(sdata[iomega, 4+4*m]), file=fi, end='')
    print('', file=fi)
fi.close()

#==============================
#     Matsubara Green function
#==============================
Gdata=np.loadtxt("Gw.dat")
nomega=Gdata.shape[0]
for iomega in range(nomega):
    m_tot=int((Gdata.shape[1]-1)/4)
    for m in range(m_tot):
      Gdata[iomega, 1+4*m] = (Gdata[iomega, 1+4*m] + Gdata[iomega, 3+4*m])/2.0
      Gdata[iomega, 3+4*m] = Gdata[iomega, 1+4*m]
      Gdata[iomega, 2+4*m] = (Gdata[iomega, 2+4*m] + Gdata[iomega, 4+4*m])/2.0
      Gdata[iomega, 4+4*m] = Gdata[iomega, 2+4*m]

fi=open("Gw.out",'w')
for iomega in range(nomega):
    m_tot=int((Gdata.shape[1]-1)/4)
    print("{:22.15f}".format(Gdata[iomega, 0]), file=fi, end='')
    for m in range(m_tot):
        print("{:20.15f}".format(Gdata[iomega, 1+4*m]), file=fi, end='')
        print("{:20.15f}".format(Gdata[iomega, 2+4*m]), file=fi, end='')
        print("{:20.15f}".format(Gdata[iomega, 3+4*m]), file=fi, end='')
        print("{:20.15f}".format(Gdata[iomega, 4+4*m]), file=fi, end='')
    print('', file=fi)
fi.close()

#==============================
# Imaginary time Green function
#==============================
Gtdata=np.loadtxt("Gt.dat")
ntau=Gtdata.shape[0]
for itau in range(ntau):
    m_tot=int((Gtdata.shape[1]-1)/2)
    for m in range(m_tot):
        Gtdata[itau, 1+2*m] = (Gtdata[itau, 1+2*m] + Gtdata[itau, 2+2*m])/2.0
        Gtdata[itau, 2+2*m] = Gtdata[itau, 1+2*m]

fi=open("Gt.out",'w')
for itau in range(ntau):
    m_tot=int((Gtdata.shape[1]-1)/2)
    print("{:22.15f}".format(Gtdata[itau, 0]), file=fi, end='')
    for m in range(m_tot):
        print("{:20.15f}".format(Gtdata[itau, 1+2*m]), file=fi, end='')
        print("{:20.15f}".format(Gtdata[itau, 2+2*m]), file=fi, end='')
    print('', file=fi)
fi.close()
