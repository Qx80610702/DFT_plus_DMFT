#!/usr/bin/env python
import numpy as np
from scipy import random
from scipy import optimize
from scipy import interpolate
from scipy import integrate
import sys, os, shutil
import maxent_routines as maxent

skrams_exe="/home/quxin/softwares/DFT_plus_DMFT/build/analy_con/skrams"

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

    Gtau = np.zeros(len(tau), dtype=float)
    df = Gm[-1].real*freq[-1]/np.pi
    print('df=', df)
    ah = FindHighFrequency(Gm,freq,Nf)
    for it,t in enumerate(tau): 
        Gtau[it] = maxent.fourpart(t,Gm,freq,ah,beta)
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

def MaximumEntropy(p, tau, Gt):
    beta = tau[-1]

    random.seed( 1 ) # seed for random numbers

    if 'x0' in p:
        omega = GiveTanMesh(p['x0'],p['L'],p['Nw'])
    else:
        omega = np.linspace(-p['L'],p['L'],2*p['Nw']+1)
    dom = np.array([0.5*(omega[1]-omega[0])]+[0.5*(omega[i+1]-omega[i-1]) for i in range(1,len(omega)-1)]+[0.5*(omega[-1]-omega[-2])])

    Gt = -Gt
    fsg=-1
    normalization = Gt[0]+Gt[-1]
    Ker = maxent.initker_fermion(omega,dom,beta,tau)
    
    print('beta=', beta)
    print('normalization=', normalization)

    # Set error
    if p['idg']:
        sxt = np.ones(len(tau))/(p['deltag']**2)
    else:
        sxt = Gt*p['deltag']
        for i in range(len(sxt)):
            if sxt[i]<1e-5: sxt[i]=1e-5
        sxt = 1./sxt**2
    
    # Set model
    if p['iflat']==0:
        model = normalization*np.ones(len(omega))/np.sum(dom)
    elif p['iflat']==1:
        model = np.exp(-omega**2/p['gwidth'])
        model *= normalization/np.dot(model,dom)
    else:
        dat=loadtxt('model.dat').transpose()
        fm=interpolate.interp1d(dat[0],dat[1])
        model = fm(omega)
        model *= normalization/np.dot(model,dom)
        #savetxt('brisi_test', vstack((tau, fsg*dot(model,Ker))).transpose())
        
    print('Model normalization=', np.dot(model,dom))

    # Set starting Aw(omega)
    Aw = random.rand(len(omega))
    Aw = Aw * (normalization/np.dot(Aw,dom))
    print('Aw normalization=', np.dot(Aw,dom))

    dlda = maxent.initdlda(omega,dom,Ker,sxt)
    
    temp=10.
    rfac=1.
    alpha=p['alpha0']
    
    for itt in range(p['Nitt']):
        print(itt, 'Restarting maxent with rfac=', rfac, 'alpha=', alpha)
        iseed = random.randint(0,sys.maxsize)
        
        maxent.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,dom,p['Asteps'],iseed)
        S = maxent.entropy(Aw,model,dom)
        Trc = maxent.lambdac(alpha,Aw,omega,dom,dlda)
        
        ratio = -2*S*alpha/Trc
        print('Finished maxent with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc, 'S=', S)
        print('   ratio=', ratio)

        np.savetxt('dos_'+str(itt), np.vstack((omega,Aw)).transpose())
        temp=0.001
        rfac=0.05
    
        if abs(ratio-1)<p['min_ratio']: break
    
        if (abs(ratio)<0.05):
            alpha *= 0.5
        else:
            alpha *= (1.+0.001*(random.rand()-0.5))/ratio
        
    for itt in range(p['Nr']):
        print('Smoothing itt ', itt)
        Aw = Broad(p['bwdth'],omega,Aw)
        Aw *= (normalization/np.dot(Aw,dom)) # Normalizing Aw
        
        np.savetxt('dos_'+str(p['Nitt']), np.vstack((omega,Aw)).transpose())
        
        temp=0.005
        rfac=0.005
        iseed = random.randint(0,sys.maxsize)
        maxent.maxent(Aw,rfac,alpha,temp,Ker,sxt,Gt,model,dom,p['Asteps'],iseed)
        
        S = maxent.entropy(Aw,model,dom)
        Trc = maxent.lambdac(alpha,Aw,omega,dom,dlda)
        ratio = -2*S*alpha/Trc
        print('Finished smoothing run with alpha=',alpha,'-2*alpha*S=',-2*alpha*S,'Trace=',Trc)
        print('   ratio=', ratio)
        
    np.savetxt('gtn', np.vstack((tau, fsg*np.dot(Aw,Ker))).transpose())
    Aw = Broad(p['bwdth'],omega,Aw)
    np.savetxt('dos.out', np.vstack((omega,Aw)).transpose())
    return (Aw, omega)

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("give input file Sigma(i\omega_n)!!!")
        exit

    # Read maxent_params.dat
    if not os.path.exists("maxent_params.dat"):
        print("Please give file maxent_prams.dat!!!")
        exit
    exec(open("maxent_params.dat").read())

    Sfile = sys.argv[1]
    Sdata = np.loadtxt(Sfile)

    Nomega = Sdata.shape[0]
    Norb = int((Sdata.shape[1]-1)/2)

    freq = Sdata[:,0]
    beta = np.pi/freq[0]
    tau = np.linspace(0, beta, params['Ntau']+1)

    for iorb in range(Norb):
        os.system("test -d orb" + str(iorb) + " && rm -rf orb" + str(iorb))
        os.mkdir("orb"+str(iorb))
        os.chdir("orb"+str(iorb))

        # Construct auxiliary function and invert Fourier transform
        Auxw = 1.0/(freq[:]*1j-Sdata[:,1+2*iorb]-Sdata[:,2+2*iorb]*1j)
        Auxt = InverseFourier(Auxw, freq, tau, beta, params['Nf'])

        # Analytic continuation of auxiliary
        (Aw, omega) = MaximumEntropy(params, tau, Auxt)
        # np.savetxt('dos', np.vstack((omega,Aw)).transpose())

        # Performs Kramars-Kronig
        izero = np.argmin(abs(omega))
        if abs(omega[izero])<1e-6:
            omega_n = np.hstack([omega[:izero],omega[izero+1:]])
            Aw_n = np.hstack([Aw[:izero],Aw[izero+1:]])
        else:
            omega_b = omega
            Aw_n = Aw
        np.savetxt('dosn', np.vstack((omega_n,Aw_n)).transpose())

        os.system(skrams_exe + " -cn 2 -s -pi dosn > Gc")
        Gcdata = np.loadtxt("Gc")
        Sfreq = Gcdata[:,0]
        Sc = Sfreq[:]-1.0/(Gcdata[:,1]+Gcdata[:,2]*1j)
        np.savetxt('sig.out', np.array([Sfreq,Sc.real,Sc.imag]).transpose())

        os.chdir("../")

    sigma0=np.loadtxt("orb0/sig.out")
    Sigma = np.zeros((sigma0.shape[0],2*Norb+1), dtype=float)
    Sigma[:,0] = sigma0[:,0]
    Sigma[:,1] = sigma0[:,1]
    Sigma[:,2] = sigma0[:,2]
    for iorb in range(1,Norb):
        sigma_tmp = np.loadtxt("orb" + str(iorb) + "/sig.out")
        Sigma[:,2*iorb+1] = sigma_tmp[:,1]
        Sigma[:,2*iorb+2] = sigma_tmp[:,2]

    sigmaf=open("Sigma.dat",'w')
    for i in range(Sigma.shape[0]):
        print("{:23.15e}".format(Sigma[i, 0]), file=sigmaf, end='')
        for iorb in range(Norb):
            print("{:25.15e}".format(Sigma[i, 2*iorb+1]), file=sigmaf, end='')
            print("{:25.15e}".format(Sigma[i, 2*iorb+2]), file=sigmaf, end='')
        print('', file=sigmaf)
    sigmaf.close()
