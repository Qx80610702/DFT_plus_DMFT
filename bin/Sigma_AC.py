#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
from scipy import integrate
import sys, os, shutil
from math import pow
from numpy.linalg import inv

maxent_exe="/home/quxin/softwares/DFT_plus_DMFT/build/maxent/maxent"

def polynomial_regression(fw, freq):
    im=complex(0.0,1.0)
    onec=complex(1.0,0.0)

    C=np.zeros((3,3),dtype=complex)
    D1=complex(0.0,0.0)
    D2=complex(0.0,0.0)
    D3=complex(0.0,0.0)

    for iomega in range(len(fw)):
        C[0,0] += complex(-1.0/pow(freq[iomega],2),0.0)
        C[0,1] += complex(0.0,1.0/pow(freq[iomega],3))
        C[0,2] += complex(1.0/pow(freq[iomega],4),0.0)
        C[1,2] += complex(0.0,-1.0/pow(freq[iomega],5))
        C[2,2] += complex(-1.0/pow(freq[iomega],6),0.0)

        D1 += fw[iomega]/(im*freq[iomega])
        D2 += -fw[iomega]/pow(freq[iomega],2)
        D3 += im*fw[iomega]/pow(freq[iomega],3)
    
    C[1,0] = C[0,1]
    C[1,1] = C[0,2]
    C[2,0] = C[0,2]
    C[2,1] = C[1,2]

    C_inver = inv(C)

    C1 = (C_inver[0,0]*D1 + C_inver[0,1]*D2 + C_inver[0,2]*D3).real
    C2 = (C_inver[1,0]*D1 + C_inver[1,1]*D2 + C_inver[1,2]*D3).real
    C3 = (C_inver[2,0]*D1 + C_inver[2,1]*D2 + C_inver[2,2]*D3).real

    return C1, C2, C3

def InverseFourier(Gm, freq, tau, beta):
    """Inverse Fourier transform which
       computes G(tau) from G(iom)
    """
    nomega = len(Gm)
    im = complex(0.0,1.0)

    C1, C2, C3 = polynomial_regression(Gm[(nomega-20):], freq[(nomega-20):])

    Gtau = np.zeros(len(tau), dtype=float)
    for it in range(len(tau)):
        t = tau[it]
        for iomega in range(nomega):
            omegan = freq[iomega]
            val_tmp = Gm[iomega] + im*C1/omegan \
                      + C2/pow(omegan,2) \
                      - im*C3/pow(omegan,3)

            Gtau[it] += ( val_tmp*(np.cos(omegan*t)
                      -im*np.sin(omegan*t)) 
                      + val_tmp.conjugate()*(np.cos(omegan*t)
                      + im*np.sin(omegan*t))
                      ).real/beta

            # val_tmp = Gm[iomega] - 1.0/(omegan*im + constC)

            # Gtau[it] += ( val_tmp*(np.cos(omegan*t)
            #           -im*np.sin(omegan*t)) 
            #           + val_tmp.conjugate()*(np.cos(omegan*t)
            #           + im*np.sin(omegan*t))
            #           ).real/beta

        Gtau[it] += -0.5*C1 + C2*(2.0*t-beta)/4.0 + C3*(beta*t-t*t)/4.0
        # Gtau[it] -= exp(constC*t)/(exp(constC*beta)+1)

    return Gtau

if __name__ == '__main__':
    #MPI initializing
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    #Parsing command line
    char_step=1
    DMFT_step=1
    if len(sys.argv) != 5:
        for step_dir in os.listdir("dmft/"):
            if step_dir.find("charge_step") != -1:
                N=int(step_dir.split("p")[-1])
                if N > char_step:
                    char_step=N
        
        for step_dir in os.listdir("dmft/charge_step" + str(char_step)):
            if step_dir.find("dmft_step") != -1:
                N=int(step_dir.split("p")[-1])
                if N > DMFT_step:
                    DMFT_step=N
    else:
        i=1
        while i < len(sys.argv):
            if sys.argv[i].lower() == "-charge_step": char_step=int(sys.argv[i+1])
            elif sys.argv[i].lower() == "-dmft_step": DMFT_step=int(sys.argv[i+1])
            else: 
                print("Unsupported parameters ", sys.argv[i])
                exit
            i = i + 2 
    
    # Number of inequivalent impuritys
    sym = open("dft/outputs_to_DMFT/symmetry.dat",'r')
    ineq_num=int(sym.readlines()[1])

    # Read maxent_params.dat
    if not os.path.exists("maxent_params.dat"):
        print("Please give file maxent_prams.dat!!!")
        exit
    exec(open("maxent_params.dat").read())

    # Chemical potential
    mu=0.0
    for lines in open("DMFT_running.log"):
        if lines.lower().find("chemical potential") != -1:
            strs=lines.strip().split(':')
            mu=float(strs[1].split('e')[0])

    # Read imaginary time self-energy
    Sw = []
    iorb2ineqm = []
    Norb=0
    for ineq in range(ineq_num):
        Sw_tmp=np.loadtxt("dmft/charge_step" + str(char_step) \
                  +"/dmft_step" + str(DMFT_step) \
                  +"/impurity" + str(ineq) + "/Sigma.dat", dtype=float)
        Norb += int((Sw_tmp.shape[1]-1)/2)
        Sw.append(Sw_tmp)
        for m in range(int((Sw_tmp.shape[1]-1)/2)):
            iorb2ineqm.append([ineq,m])

    #Construct folder
    if rank ==0:
      if os.path.exists("self-energy"):
          shutil.rmtree("self-energy")
      os.mkdir("self-energy")
      for imp in range(ineq_num):
          os.mkdir("self-energy/impurity"+str(imp))
    comm.Barrier()


    for iorb in range(Norb):
        if iorb%size!=rank :
            continue

        imp = iorb2ineqm[iorb][0]
        mag = iorb2ineqm[iorb][1]

        os.mkdir("self-energy/impurity"+str(imp)+"/orb"+str(mag))
        os.chdir("self-energy/impurity"+str(imp)+"/orb"+str(mag))

        Nomega = Sw[imp].shape[0]
        freq = Sw[imp][:,0]
        beta = np.pi/freq[0]
        tau = np.linspace(0, beta, params['Ntau']+1)
        
        #====================Analytical continuation======================
        #   Construct auxiliary function and invert Fourier transform
        # ================================================================
        Sigoo = Sw[imp][-1,1+2*mag]+Sw[imp][-1,2+2*mag]*1.0j
        # Sigoo = Sw[imp][-1,1+2*mag]
        Auxw = 1.0/(freq[:]*1.0j-Sw[imp][:,1+2*mag]-Sw[imp][:,2+2*mag]*1.0j + Sigoo + mu)
        Auxt = InverseFourier(Auxw, freq, tau, beta)

        # #write Gtau.dat
        Gf=open("Gtau.dat",'w')
        for i in range(tau.shape[0]):
            print("{:23.15e}".format(tau[i]), file=Gf, end='')
            print("{:25.15e}".format(Auxt[i]), file=Gf)
        Gf.close()

        #write maxent_input.dat
        ofs=open("maxent_input.dat", 'w')
        ofs.write("%d         #Nt:Number of time points\n"%tau.shape[0])
        ofs.write("%d         #idg:error scheme: idg=1 -> sigma=deltag ; idg=0 -> sigma=deltag*G(tau)\n"%params['idg'])
        ofs.write(str(params['deltag']) + "         #deltag:error\n")
        ofs.write("%.15f          #beta:inverse temperature\n"%tau[-1])
        ofs.write("%d         #Nw:number of frequencies\n"%params['Nw'])
        ofs.write(str(params['Dw']) + "         #Dw:delta of frequency\n")
        ofs.write("%d         #Asteps:number of annealig steps\n"%params['Asteps'])
        ofs.write(str(params['alpha0']) + "         #alpha0:starting alpha\n")
        ofs.write("%d         #Nr:number of smooth runs\n"%params['Nr'])
        ofs.write("%d         #iflat:iflat=0 : constant model, iflat=1 : gaussian of width gwidth, iflat=2 : input using file model.dat\n"%params['iflat'])
        ofs.close()

        #Analytical continuation
        os.system(maxent_exe + " 1>running.log 2>running.error")

        #=======================Performs Kramars-Kronig transform===================
        ieta = float(params['eta'])*1j
        Aux_im = np.loadtxt("dos")
        Aux_omega = np.zeros(Aux_im.shape[0], dtype=complex)
        for i in range(Aux_omega.shape[0]):
            omega = Aux_im[i,0]
            ker = Aux_im[:,1]/(omega - Aux_im[:,0] + ieta)
            Aux_omega[i]=integrate.trapz(ker[:], Aux_im[:,0] )
        
        Sc = Aux_im[:,0] - 1.0/Aux_omega[:] + Sigoo + mu
        np.savetxt('sig.out', np.array([Aux_im[:,0], Sc.real, Sc.imag]).transpose())

        os.chdir("../../../")

    comm.Barrier()

    if rank==0 :
        for imp in range(ineq_num):
            os.chdir("self-energy/impurity"+str(imp))

            sigma0=np.loadtxt("orb0/sig.out")
            mag_num=0
            for dir in os.listdir("./"):
                if dir.find("orb") != -1:
                    mag_num = mag_num + 1
            Sigma_data = np.zeros((sigma0.shape[0],2*mag_num+1), dtype=float)
            Sigma_data[:,0] = sigma0[:,0]
            Sigma_data[:,1] = sigma0[:,1]
            Sigma_data[:,2] = sigma0[:,2]
            for m in range(1,mag_num):
                sigma_tmp = np.loadtxt("orb" + str(m) + "/sig.out")
                Sigma_data[:,2*m+1] = sigma_tmp[:,1]
                Sigma_data[:,2*m+2] = sigma_tmp[:,2]

            sigmaf=open("Sigma_omega.dat",'w')
            for i in range(Sigma_data.shape[0]):
                print("{:23.15e}".format(Sigma_data[i, 0]), file=sigmaf, end='')
                for m in range(mag_num):
                    print("{:25.15e}".format(Sigma_data[i, 2*m+1]), file=sigmaf, end='')
                    print("{:25.15e}".format(Sigma_data[i, 2*m+2]), file=sigmaf, end='')
                print('', file=sigmaf)
            sigmaf.close()

            os.chdir("../../")
