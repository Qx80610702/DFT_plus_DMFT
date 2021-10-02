#!/usr/bin/env python
import numpy as np
import os, sys
from Sigma_AC import MaximumEntropy

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("give input imeginary time Green's function!!!")
        exit
    
    # Read maxent_params.dat
    if not os.path.exists("maxent_params.dat"):
        print("Please give file maxent_prams.dat!!!")
        exit
    exec(open("maxent_params.dat").read())

    Gt_data=np.loadtxt(sys.argv[1])
    Ntau=Gt_data.shape[0]
    Norb=Gt_data.shape[1]-1

    for iorb in range(Norb):
        os.system("test -d orb" + str(iorb) + " && rm -rf orb" + str(iorb))
        os.mkdir("orb"+str(iorb))
        os.chdir("orb"+str(iorb))

        # Analytic continuation of auxiliary
        (Aw, omega) = MaximumEntropy(params, Gt_data[:,0], Gt_data[:,iorb+1])

        os.chdir("../")

    dos0=np.loadtxt("orb0/dos.out")
    Aw_data = np.zeros((dos0.shape[0],Norb+1), dtype=float)
    for iorb in range(Norb):
        dos_tmp=np.loadtxt("orb" + str(iorb) + "/dos.out")
        if(iorb==0):
            for i in range(dos_tmp.shape[0]):
                Aw_data[i,0] = dos_tmp[i,0]
                Aw_data[i,iorb+1] = dos_tmp[i,1]
        else:
            for i in range(dos_tmp.shape[0]):
                Aw_data[i,iorb+1] = dos_tmp[i,1]

    Awf=open("spectrum.dat",'w')
    for i in range(Aw_data.shape[0]):
        print("{:23.15e}".format(Aw_data[i, 0]), file=Awf, end='')
        for iorb in range(Norb):
            print("{:25.15e}".format(Aw_data[i, iorb+1]), file=Awf, end='')
        print('', file=Awf)
    Awf.close()


