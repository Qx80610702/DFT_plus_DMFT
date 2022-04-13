#!/usr/bin/env python
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import ceil
import sys

def read_latt_vector(DFT_solver):
    latvec = []
    if DFT_solver.lower() == "aims":
        for line in open("./dft/geometry.in"):
            line = line.strip().split("#")[0]
            words = line.strip().split()
            if len(words) == 0:
                continue
            if words[0] == "lattice_vector":
                if len(words) != 4:
                    print("geometry.in: Syntax error in line '"+line+"'")
                    exit()
                latvec += [ np.array(list(map(float,words[1:4]))) ]
    elif DFT_solver.lower() == "abacus":
        latt=1.0
        STRU=open("./dft/STRU",'r')
        lines=STRU.readlines()
        i=0
        while i < len(lines):
            words = lines[i].strip().split()
            if len(words) == 0:
                i += 1
                continue
            if words[0].lower() == "lattice_constant":
                latt=float(lines[i+1].split()[0])
                i += 2
                continue
            elif words[0].lower() == "lattice_vectors":
                latvec += [ latt*np.array(list(map(float,lines[i+1].strip().split()))) ]
                latvec += [ latt*np.array(list(map(float,lines[i+2].strip().split()))) ]
                latvec += [ latt*np.array(list(map(float,lines[i+3].strip().split()))) ]
                break
            else:
                i += 1
    else:
        print("Unsuported DFT solver!!!")
        exit()
    
    return np.asarray(latvec)

def cal_rec_latvec(latvec):
    rlatvec = []
    volume = (np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])))
    rlatvec.append(np.array(2*np.pi*np.cross(latvec[1,:],latvec[2,:])/volume))
    rlatvec.append(np.array(2*np.pi*np.cross(latvec[2,:],latvec[0,:])/volume))
    rlatvec.append(np.array(2*np.pi*np.cross(latvec[0,:],latvec[1,:])/volume))

    return np.asarray(rlatvec)

def parsing_DFT_solver():
    DFT_solver = "aims"  #default value
    for line in open("DMFT.in",'r'):
        words = line.strip().split()
        if len(words) == 0:
            continue
        elif words[0].lower() == "dft_solver":
            DFT_solver = words[1]

    return DFT_solver.lower()

def read_kpath():
    kpath = []
    knames = []
    nks = 0
    for line in open("kpath.in",'r'):
        words = line.strip().split()
        if len(words)==0:
            continue
        if words[0].lower()=="path":
            if len(words) != 5 :
              print("kpath.in: Syntax error in line '"+line+"'")
              exit
            kpath.append(list(words[1:4]))
            knames.append(words[4])
        elif words[0].lower()=="total_kpoints":
            nks=int(words[1])

    return kpath, knames, nks

def cal_kpoints(kpath, knames, nks, rlatvec):
    kpoints = []
    sym_kpt = {}
    band_totlength = 0.0 # total length of all band segments
    for i in range(len(kpath)-1):
        start = np.array(list(map(float,kpath[i])))
        end = np.array(list(map(float,kpath[i+1])))
        length = norm(np.dot(rlatvec,end) - np.dot(rlatvec,start))
        band_totlength += length
    
    delta = band_totlength/nks
    
    count=0
    for i in range(len(kpath)-1):
        start = np.array(list(map(float,kpath[i])))
        end = np.array(list(map(float,kpath[i+1])))
        length = norm(np.dot(rlatvec,end) - np.dot(rlatvec,start))
        nks_seg = ceil(length/delta)
        seg_kx = np.linspace(start[0], end[0], nks_seg, dtype=float)
        seg_ky = np.linspace(start[1], end[1], nks_seg, dtype=float)
        seg_kz = np.linspace(start[2], end[2], nks_seg, dtype=float)
        if i==0:
          sym_kpt['0'] = knames[0]
          for ik_seg in range(nks_seg):
              kpoints.append([seg_kx[ik_seg],seg_ky[ik_seg],seg_kz[ik_seg]])
          count += nks_seg
          sym_kpt[str(count-1)] = knames[1]
        else:
          for ik_seg in range(1,nks_seg):
              kpoints.append([seg_kx[ik_seg],seg_ky[ik_seg],seg_kz[ik_seg]])
          count += nks_seg-1
          sym_kpt[str(count-1)] = knames[i+1]

    return kpoints, sym_kpt
        
if __name__ == '__main__':
    matplotlib.rcParams['lines.linewidth'] = 0.5
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'out'
    # matplotlib.rcParams['figure.dpi'] = 500
    # matplotlib.rcParams['savefig.dpi'] = 500

    #=======parameters========
    Awf = " "
    eup = 20.0
    edw = -20.0
    Etciks = 2.0
    intensity = 1.0
    showcolorbar = 1
    outfig = "bands.png"

    #======frequency====
    omega = np.loadtxt("self-energy/impurity0/Sigma_omega.dat")[:,0]
    eup=omega[-1]
    edw=omega[0]

    # ====pasrsing command line======
    if(len(sys.argv)<2):
        print("Please give spectrum file")
        exit
    else:
        i=1
        while i < len(sys.argv):
            if sys.argv[i].find('-') == -1:
                Awf=sys.argv[i]
                i = i + 1
            else:
                if sys.argv[i].lower() == "-eup": eup=float(sys.argv[i+1])
                elif sys.argv[i].lower() == "-edw": edw=float(sys.argv[i+1])
                elif sys.argv[i].lower() == "-eticks": Etciks=float(sys.argv[i+1])
                elif sys.argv[i].lower() == "-contrast": intensity=float(sys.argv[i+1])
                elif sys.argv[i].lower() == "-colorbar": showcolorbar=int(sys.argv[i+1])
                elif sys.argv[i].lower() == "-o": outfig=sys.argv[i+1]
                else: 
                    print("Unsupported parameters ", sys.argv[i])
                    exit
                i = i + 2
    
    if Awf == " " :
        print("Please give spectrum file")
        exit

    Awk = np.loadtxt(Awf)

    max_index = 0
    for i in range(len(omega)):
        if omega[i] >= eup:
          max_index=i
          break
    
    min_index = len(omega)
    for i in range(len(omega)):
        if omega[i] >= edw:
          min_index=i
          break

    omegaw = omega[min_index:max_index+1]
    Awkw = Awk[min_index:max_index+1,:]

    #======Parsing high symmetry k-points======
    DFT_solver = parsing_DFT_solver()
    latvec = read_latt_vector(DFT_solver)
    rlatvec = cal_rec_latvec(latvec)
    kpath, knames, nks = read_kpath()
    kpoints, sym_kpt = cal_kpoints(kpath, knames, nks, rlatvec)

    vmm = [0,max(map(max,Awkw))*intensity]

    plt.xticks([])
    plt.imshow(Awkw, origin='lower', interpolation='bilinear', cmap=cm.plasma, vmin=vmm[0], vmax=vmm[1], aspect='auto')
    if showcolorbar==1: plt.colorbar()

    #===== x axis=====
    tickx = []
    tickn = []
    for ik, names in sym_kpt.items():
        tickx += [ int(ik) ]
        if names.lower()=="gamma":
            names = "$\\Gamma$"
        elif names.lower()=="delta":
            names = "$\\Delta$"
        elif names.lower()=="lambda":
            names = "$\\Lambda$"
        elif names.lower()=="sigma":
            names = "$\\Sigma$"
        tickn += [ names ]
        plt.axvline(float(ik), color='w',linestyle="-")
    plt.xticks(tickx, tickn)

    #====y axis====
    exec(open("maxent_params.dat").read())
    npoints=int(Etciks/float(params['Dw']))
    plt.yticks(np.arange(0,len(omegaw),npoints), omegaw[::npoints])

    if eup>0.0 and edw <0.0:
        Fermi = int(abs(edw)/float(params['Dw']))
        plt.axhline(Fermi, color='w', linestyle=":")

    plt.savefig(outfig)
    plt.show()
