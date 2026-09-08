from sys import argv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

if (len(argv) == 5):
    prefix = argv[1]
    ntime = int(argv[2])
    nz = int(argv[3])
    ny = int(argv[4])

_files = [(4, prefix+'_interface.S'),
          (6, prefix+'_interface.Vhat')]

for nfields, fname in _files:

    sdat = np.empty((ntime, nfields, nz, ny), dtype='float32')

    f = open('data/'+fname, 'rb')
    temp_dat = np.fromfile(f, 'float32')
    sdat = np.reshape(temp_dat[:sdat.size], sdat.shape)

    for i in range(0, ntime, ntime/10):
        fig = plt.figure(figsize=(8.0, (8.0*nz)/ny))
        plt.suptitle(fname + ' at ' + str(i))

        ax1 = fig.add_subplot(111)
        z = sdat[i, 2, :, :]
        levels = MaxNLocator(nbins=51).bin_boundaries(z.min(), min(10.0,z.max()))
        cmap = plt.get_cmap('RdBu')
        norm = BoundaryNorm(levels, ncolors=cmap.N)
        plt.pcolormesh(sdat[i, 2, :, :],cmap=cmap, norm=norm)
        plt.axis([0, ny-1, 0, nz-1])
        #plt.colorbar()
        ax1.set_aspect('auto')

        plt.savefig(fname + '_{0:08d}.png'.format(i))
        plt.close()
