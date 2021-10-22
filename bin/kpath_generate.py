#!/usr/bin/env python

import numpy as np
from numpy.linalg import norm
from math import ceil

def read_latt_vector(DFT_solver):
    latvec = []
    if DFT_solver.lower() == "aims":
        for line in open("../DFT/geometry.in"):
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
        print("abacus")
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

def out_kpoints_to_dft_solver(DFT_solver,kponits):
    kweight = float(1.0/len(kpoints))
    if DFT_solver == "aims":
        ofs = open("../DFT/k_list.in", 'w')
        print("1 1 1", file=ofs)
        print(len(kpoints), file=ofs)
        for ik in range(len(kpoints)):
            print("{:20.15f}".format(kpoints[ik][0]), file=ofs, end='')
            print("{:20.15f}".format(kpoints[ik][1]), file=ofs, end='')
            print("{:20.15f}".format(kpoints[ik][2]), file=ofs, end='')
            print("{:20.15f}".format(kweight), file=ofs)
        ofs.close()

if __name__ == '__main__':
    DFT_solver = parsing_DFT_solver()
    latvec = read_latt_vector(DFT_solver)
    rlatvec = cal_rec_latvec(latvec)
    kpath, knames, nks = read_kpath()
    kpoints, sym_kpt = cal_kpoints(kpath, knames, nks, rlatvec)
    out_kpoints_to_dft_solver(DFT_solver, kpoints)





