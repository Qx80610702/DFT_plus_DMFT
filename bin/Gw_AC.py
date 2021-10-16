#!/usr/bin/env python
import numpy as np
import os, sys, shutil
from mpi4py import MPI
from Sigma_AC import Broad

maxent_exe="/home/quxin/softwares/DFT_plus_DMFT/build/analy_continuation/maxent"

if __name__ == '__main__':
    #MPI initializing
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    #Parsing command line
    if len(sys.argv)<2:
        DMFT_step=1
        for step_dir in os.listdir("impurity_solving/"):
            N=int(step_dir[-1])
            if N > DMFT_step:
                DMFT_step=N
    else: 
        DMFT_step = int(sys.argv[1])
    
    #Number of inequivalent impuritys
    sym = open("../DFT/outputs_to_DMFT/symmetry.dat",'r')
    ineq_num=int(sym.readlines()[1])

    # Read maxent_params.dat
    if not os.path.exists("maxent_params.dat"):
        print("Please give file maxent_prams.dat!!!")
        exit
    exec(open("maxent_params.dat").read())

    # Read imaginary time Green's function 
    Gt_data = []
    iorb2ineqm = []
    Norb=0
    for ineq in range(ineq_num):
        Gt_tmp=np.loadtxt("impurity_solving/step" \
                        +str(DMFT_step)+"/impurity"\
                        +str(ineq)+"/Gt.dat", dtype=float)
        Norb += Gt_tmp.shape[1]-1
        Gt_data.append(Gt_tmp)
        for m in range(Gt_tmp.shape[1]-1):
            iorb2ineqm.append([ineq,m])

    #Construct folder
    if rank ==0:
      if os.path.exists("Aw_loc"):
          shutil.rmtree("Aw_loc")
      os.mkdir("Aw_loc")
      for imp in range(ineq_num):
          os.mkdir("Aw_loc/impurity"+str(imp))
    comm.Barrier()

    for iorb in range(Norb):
        if iorb%size!=rank : 
            continue

        imp = iorb2ineqm[iorb][0]
        mag = iorb2ineqm[iorb][1]

        os.mkdir("Aw_loc/impurity"+str(imp)+"/orb"+str(mag))
        os.chdir("Aw_loc/impurity"+str(imp)+"/orb"+str(mag))

        Ntau=Gt_data[imp].shape[0]
        Ntau_tmp = Ntau
        if Ntau>800 and Ntau%2==1:
            Ntau_tmp = int(Ntau/2) + 1
            Gtau = np.zeros((Ntau_tmp,2), dtype=float)
            for itau in range(Ntau_tmp):
                Gtau[itau,0] = Gt_data[imp][2*itau,0]
                Gtau[itau,1] = Gt_data[imp][2*itau,mag+1]
        else :
            Gtau = np.zeros((Ntau_tmp,2),dtype=float)
            for itau in range(Ntau_tmp):
                Gtau[itau,0] = Gt_data[imp][itau,0]
                Gtau[itau,1] = Gt_data[imp][itau,mag+1]

        #write maxent_input.dat
        ofs=open("maxent_input.dat", 'w')
        ofs.write("%d         #Nt:Number of time points\n"%Gtau.shape[0])
        ofs.write("%d         #idg:error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)\n"%params['idg'])
        ofs.write(str(params['deltag']) + "         #deltag:error\n")
        ofs.write("%.15f          #beta:inverse temperature\n"%Gtau[-1,0])
        ofs.write("%d         #Nw:number of frequencies\n"%params['Nw'])
        ofs.write(str(params['Dw']) + "         #Dw:delta of frequency\n")
        ofs.write("%d         #Asteps:number of annealig steps\n"%params['Asteps'])
        ofs.write(str(params['alpha0']) + "         #alpha0:starting alpha\n")
        ofs.write("%d         #Nr:number of smooth runs\n"%params['Nr'])
        ofs.write("%d         #iflat:iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat\n"%params['iflat'])
        ofs.close()

        #write Gtau.dat
        Gf=open("Gtau.dat",'w')
        for i in range(Gtau.shape[0]):
            print("{:23.15e}".format(Gtau[i, 0]), file=Gf, end='')
            print("{:25.15e}".format(Gtau[i, 1]), file=Gf)
        Gf.close()

        #Analytical continuation
        os.system(maxent_exe + " 1>running.log 2>running.error")
        #Smothing
        DOS=np.loadtxt("dos")
        Aw = Broad(params['bwdth'], DOS[:,0], DOS[:,1])

        dosf=open("dos.out",'w')
        for i in range(DOS.shape[0]):
            print("{:23.15e}".format(DOS[i, 0]), file=dosf, end=' ')
            print("{:25.15e}".format(Aw[i]), file=dosf)
        dosf.close()

        os.chdir("../../../")

    comm.Barrier()

    if rank==0 :
        for imp in range(ineq_num):
            os.chdir("Aw_loc/impurity"+str(imp))

            dos0=np.loadtxt("orb0/dos.out")
            mag_num = len(os.listdir("./"))
            Aw_data = np.zeros((dos0.shape[0],mag_num+1), dtype=float)
            for m in range(mag_num):
                dos_tmp=np.loadtxt("orb" + str(m) + "/dos.out")
                if(m==0):
                    for i in range(dos_tmp.shape[0]):
                        Aw_data[i,0] = dos_tmp[i,0]
                        Aw_data[i,m+1] = dos_tmp[i,1]
                else:
                    for i in range(dos_tmp.shape[0]):
                        Aw_data[i,m+1] = dos_tmp[i,1]

            Awf=open("spectrum.dat",'w')
            for i in range(Aw_data.shape[0]):
                print("{:23.15e}".format(Aw_data[i, 0]), file=Awf, end='')
                for m in range(Aw_data.shape[1]-1):
                    print("{:25.15e}".format(Aw_data[i, m+1]), file=Awf, end='')
                print('', file=Awf)
            Awf.close()

            os.chdir("../../")

