from sys import argv
import numpy as np
import matplotlib.pyplot as plt


dir1 = 'truth'
dir2 = 'data'
prefix = 'test_rup_curv_fric'
ngrid = 21
ntime = 46

if (len(argv) ==6):
    dir1 = argv[1]
    dir2 = argv[2]
    prefix = argv[3]
    ngrid = int(argv[4])
    ntime = int(argv[5])

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

    filenames = [dir1+'/'+fname, dir2+'/'+fname]

    for n, name in enumerate(filenames):
        f = open(name, 'rb')
        temp_dat = np.fromfile(f, 'float32')
        fields[n] = np.reshape(temp_dat, fields.shape[1:])

    diff = fields - fields[0]

    for i in [0, ntime-1]:
        fig = plt.figure(figsize=(15, 4))
        plt.suptitle(fname + ' at ' + str(i))

        ax1 = fig.add_subplot(131)
        plt.pcolormesh(diff[1, i, 2, :, :])
        plt.axis([0, ngrid, 0, ngrid])
        plt.colorbar()
        plt.title('Difference')
        ax1.set_aspect('equal')

        ax2 = fig.add_subplot(132)
        plt.pcolormesh(fields[0, i, 2, :, :])
        plt.axis([0, ngrid, 0, ngrid])
        plt.colorbar()
        plt.title('Truth')
        ax2.set_aspect('equal')

        ax4 = fig.add_subplot(133)
        plt.pcolormesh(fields[1, i, 2, :, :])
        plt.axis([0, ngrid, 0, ngrid])
        plt.colorbar()
        plt.title('Test')
        ax4.set_aspect('equal')

        plt.savefig(fname + '.' + str(i) + '.png')
        # plt.show()
        plt.close()
