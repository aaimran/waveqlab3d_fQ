#!/usr/bin/env python
from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt


dir1 = 'truth'
dir2 = 'data'
ngrid = 21
ntime = 46

prefix = argv[2]
path = argv[1]

_files = [(4, prefix+'_interface.S'),
          (3, prefix+'_interface.Svel'),
          # (1, prefix+'_interface.trup'),
          # (1, prefix+'_interface.state'),
          (9, prefix+'_interface.Uface'),
          (9, prefix+'_interface.Vface'),
          (6, prefix+'_interface.Uhat'),
          (6, prefix+'_interface.Vhat')]

for nfields, fname in _files:
    fields = np.empty((2, ntime, nfields, ngrid, ngrid), dtype='float32')

    filenames = [path+'/'+dir1+'/'+fname, path+'/'+dir2+'/'+fname]

    for n, name in enumerate(filenames):
        f = open(name, 'rb')
        temp_dat = np.fromfile(f, 'float32')
        fields[n] = np.reshape(temp_dat, fields.shape[1:])

    diff = fields - fields[0]
    diffabs = np.sqrt(diff*diff)

    if np.amax(diffabs) > 1.0e-3:
      print fname, np.amax(diffabs)
      exit(1)