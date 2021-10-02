#!/usr/bin/env python

import h5py
import math

def read_param(h5, name):
  if '/parameters/dictionary/'+name in h5:
    return h5['/parameters/dictionary/'+name].value
  elif '/parameters/'+name in h5:
    return h5['/parameters/'+name].value
  else:
    raise RuntimeError("Parameter "+ name + " not found")

def read_sp_gf(h5, path, sites, spins, out_file):
  N = h5[path].shape[0]
  M = h5[path].shape[1]
  if M != sites*spins:
    raise RuntimeError("Error in reading the single particle Green function from the output hdf5 file of QMC solver")
  data = h5[path].value  #A multiple list
  #print(data.shape)

  f=open(out_file,'w')
  #f.write(str(N) + "  " + str(sites) + "  " + str(spins) + "\n")
  for n in range(N):
    for i in range(M):
      for j in range(M):
        re = '{:.15f}'.format(data[n][i][j][0])
        im = '{:.15f}'.format(data[n][i][j][1])

        f.write(
          str(n) + "  "
          + str(i) + "  "
          + str(j) + "  "
          + str(re) + "  "
          + str(im) + "\n"
        )

  f.close()

def read_h5(p):
  r = {}
  print('Reading '+ p +'.out.h5')
  h5 = h5py.File(p+'.out.h5','r')

  r["SITES"] = read_param(h5, 'model.sites')
  r['spins'] = read_param(h5, 'model.spins')
  r["BETA"] = read_param(h5, 'model.beta')

  r["Gtau"] = read_sp_gf(h5,'gtau/data', r["SITES"], r['spins'], 'G_tau.dat')
  r["Gomega"] = read_sp_gf(h5,'gf/data', r["SITES"], r['spins'], 'G_omega.dat')

if __name__ == "__main__":
  read_h5("input")
