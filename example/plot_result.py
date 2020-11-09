import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# ***************************************************

class data:
    """
    unpacks data into named variables
    """
    def __init__(self, cub):
        self.z = cub[:,:,0].squeeze()
        self.ltau = cub[:,:,1].squeeze()
        self.Tg = cub[:,:,2].squeeze()
        self.Vz = cub[:,:,3].squeeze()
        self.Vt = cub[:,:,4].squeeze()
        self.Bx = cub[:,:,5].squeeze()
        self.By = cub[:,:,6].squeeze()
        self.Bz = cub[:,:,7].squeeze()
        self.Ne = cub[:,:,8].squeeze()
        self.rho = cub[:,:,0].squeeze()
        self.nH = cub[:,:,10::].squeeze()

# ***************************************************


if __name__ == "__main__":
    plt.close("all"); plt.ion()
    
    #
    # Read data
    #
    dorig = data(fits.open("FALC_5x5.fits",'readonly')[0].data[:]*1)
    dopti = data(fits.open("FALC_5x5_optimized.fits",'readonly')[0].data[:]*1)


    
    #
    # make plot
    #
    f, ax = plt.subplots(nrows=1, ncols=1, figsize=(7,3))
    
    ax.plot(dorig.z[0,0]*1.e-5, dorig.Tg[0,0]*1.e-3, '.-', color='black', label='original')
    ax.plot(dopti.z[0,0]*1.e-5, dopti.Tg[0,0]*1.e-3, '+', color='orangered', label='depth-opt')

    ax.set_ylim(4, 15)
    ax.set_ylabel('T [kK]')
    ax.set_xlabel('z [km]')
    
    ax.legend(frameon = True, loc = 'upper left')
    
    f.set_tight_layout(True)
    f.show()
    
