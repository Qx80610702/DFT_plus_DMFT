#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
from scipy import random
from scipy import optimize
from scipy import interpolate
from scipy import integrate
import sys, os, shutil

skrams_exe="/home/quxin/softwares/DFT_plus_DMFT/build/analy_continuation/skrams"
maxent_exe="/home/quxin/softwares/DFT_plus_DMFT/build/analy_continuation/maxent"

def InverseFourier(Gm, freq, tau, beta, Nf=40):
    """Inverse Fourier transform which
       computes G(tau) from G(iom)
    """
    def FindHighFrequency(Gm,om,Nf):
        S=0.; Sx=0.; Sy=0.; Sxx=0.; Sxy=0;
        for j in range(len(om)-Nf,len(om)):
            x = om[j]
            y = Gm[j].imag * x
            x2= x**2
            Sy += y
            Sx += 1/x2
            Sxx += 1/x2**2
            Sxy += y/x2
            S += 1
    
        dd = S*Sxx-Sx**2
        a = (Sxx*Sy-Sx*Sxy)/dd
        bx = (S*Sxy-Sx*Sy)/dd
        ah = -a;
        if abs(ah-1.0)<1e-3: ah=1.0
        return ah

    def fourpart(t,Gm,freq,ah,beta):
        sum=0.0
        Nw = freq.shape[0]
        for iw in range(Nw):
            sum += np.cos(freq[iw]*t)*Gm[iw].real + \
                   np.sin(freq[iw]*t)*( Gm[iw].imag+ah/freq[iw] )
        return 2.0*sum/beta-0.5*ah

    Gtau = np.zeros(len(tau), dtype=float)
    df = Gm[-1].real*freq[-1]/np.pi
    print('df=', df)
    ah = FindHighFrequency(Gm,freq,Nf)
    for it,t in enumerate(tau): 
        Gtau[it] = fourpart(t,Gm,freq,ah,beta)
    Gtau[0] += df
    Gtau[-1] -= df

    return Gtau

def GiveTanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return np.array([L-w/np.tan(d), x0-w*np.tan(np.pi/(2*Nw)-d/Nw) ])
    
    xi=x0/L
    d0 = Nw/2.*(np.tan(np.pi/(2*Nw))-np.sqrt(np.tan(np.pi/(2*Nw))**2 - 4*xi/Nw))
    w0 = L*d0

    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    om = w*np.tan(np.linspace(0,1,2*Nw+1)*(np.pi-2*d) -np.pi/2+d)
    return om

def Broad(width, om, fw):
    " Broadens the data with gaussian of width=width"
    def MakeTanMesh(N, tanc, tanw, b0, b1):
        if not(b0<b1): print("Relation must hold: b0<b1!")
        if not(b0<tanw and tanw<b1): print("Relation mesu hold: b0<tanw<b1!")
        if not(b0>0): print("b0 must be positive!")
        du = np.arctan(((tanc-b0)/tanw))
        b1n = np.arctan((b1-tanc)/tanw)+du
        m0 = [tanc + tanw * np.tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
        return np.hstack( (-np.array(m0[::-1]), np.array([0]+m0) ) )

    if width<1e-5: return fw
    
    w=width
    x = MakeTanMesh(200,0.0,w,w/50,w*20)
    fwi = interpolate.interp1d(om, fw)
    fwn=[]
    for im in range(len(om)):
        x1 = list(filter(lambda t: t>=om[im]-x[-1] and t<=om[im]-x[0], om))
        x2 = list(filter(lambda t: t>=om[0] and t<=om[-1], x+om[im]))
        eps = sorted(np.hstack((x1, x2)))
        x3 = om[im]-eps
        gs = np.exp(-x3**2/(2*w**2))/(np.sqrt(2*np.pi)*w)
        yn = integrate.trapz(fwi(eps) * gs, x=eps)
        fwn.append(yn)
    return np.array(fwn)

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

    # Read imaginary time self-energy
    Sw = []
    iorb2ineqm = []
    Norb=0
    for ineq in range(ineq_num):
        Sw_tmp=np.loadtxt("impurity_solving/step" \
                        +str(DMFT_step)+"/impurity"\
                        +str(ineq)+"/Sigma.dat", dtype=float)
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
        
        #=============Analytical continuation=============================
        # Construct auxiliary function and invert Fourier transform
        Auxw = 1.0/(freq[:]*1j-Sw[imp][:,1+2*mag]-Sw[imp][:,2+2*mag]*1j)
        Auxt = InverseFourier(Auxw, freq, tau, beta, params['Nf'])

        #write Gtau.dat
        Gf=open("Gtau.dat",'w')
        for i in range(tau.shape[0]):
            print("{:23.15e}".format(tau[i]), file=Gf, end='')
            print("{:25.15e}".format(Auxt[i]), file=Gf)
        Gf.close()

        Gtau = np.loadtxt("Gtau.dat")

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

        #=======================Performs Kramars-Kronig===================
        Aux_omega = np.loadtxt("dos.out")
        omega=Aux_omega[:,0]
        izero = np.argmin(abs(omega))
        if abs(omega[izero])<1e-6:
            omega_n = np.hstack([omega[:izero],omega[izero+1:]])
            Aw_n = np.hstack([Aux_omega[:izero,1],Aux_omega[izero+1:,1]])
        else:
            omega_b = omega
            Aw_n = Aux_omega[:,1]
        np.savetxt('dosn', np.vstack((omega_n,Aw_n)).transpose())

        os.system(skrams_exe + " -cn 2 -s -pi dosn > Gc")
        Gcdata = np.loadtxt("Gc")
        Sfreq = Gcdata[:,0]
        Sc = Sfreq[:]-1.0/(Gcdata[:,1]+Gcdata[:,2]*1j)
        np.savetxt('sig.out', np.array([Sfreq,Sc.real,Sc.imag]).transpose())

        os.chdir("../../../")

    comm.Barrier()

    if rank==0 :
        for imp in range(ineq_num):
            os.chdir("self-energy/impurity"+str(imp))

            sigma0=np.loadtxt("orb0/sig.out")
            mag_num = len(os.listdir("./"))
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
