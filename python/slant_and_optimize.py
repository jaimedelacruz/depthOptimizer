import simulation_tools as st
from importlib import reload
import numpy as np
from astropy.io import fits
import os
import helita.sim.rh15d as rh

reload(st)

M_TO_CM = 100.0
M3_TO_CM3 = M_TO_CM**3
KG_TO_GR  = 1000.
RHO_TRANS = KG_TO_GR / M3_TO_CM3
TESLA_TO_GAUSS = 10000.

#
# Some wrappers to write/read FITS files
#

# ********************************************************* 

def writefits(name, d):
    io = fits.PrimaryHDU(d)
    io.writeto(name, overwrite=True)
    
# ********************************************************* 
    
def readfits(name, ext=0):
    return fits.getdata(name, ext=ext)

# ********************************************************* 

if __name__ == "__main__":
    #
    # Let's assume that we have the following variables 
    # with dimensions (nz, ny, nx), in python axis ordering.
    # In IDL the equivalent would be (nx, ny, nz). All in CGS units.
    # Tg, Vz, Vx, Vy, Bz, By, Bx, rho, Ne and potentially nH(ny, nx, 6, nz).
    # If we don't have nH we can leave it blank and the code will fill it.
    #
    # Additionally we have z, x and y, which are the grid sampling in the z, x, y direction (1D arrays).
    #
    # 


    #
    # If we want to slant and project the model to a different mu angle
    #
    mu_x = 1.0
    mu_y = 0.65

    if(np.abs(mu_x) < 0.999 or np.abs(mu_y) < 0.999):
        #
        # Init slanter object
        #
        slanter = st.SlantModel(z, y, x, mu_y = mu_y, mu_x = mu_x)

        #
        # slant all variables
        #
        Tg  = slanter.slantVariable(Tg)
        Vz  = slanter.slantVariable(Vz)
        Vy  = slanter.slantVariable(Vy)
        Vx  = slanter.slantVariable(Vx)
        Bz  = slanter.slantVariable(Bz)
        By  = slanter.slantVariable(By)
        Bx  = slanter.slantVariable(Bx)
        rho = slanter.slantVariable(rho)
        Ne  = slanter.slantVariable(Ne)

        #
        # Get new z-scale and copy it to a 3D cube
        #
        znew = slanter.get_z_nrew()

        z3d = np.zeros(Tg.size, dtype='float32')
        for ii in range(z.size):
            z3d[ii,:,:] = znew[ii]


        #
        # project B-field and V-field into the new L.O.S
        #
        slanter.project_field(Vy, Vx, Vz)
        slanter.project_field(By, Bx, Bz)
    else: # if no slanting is required, just create a 3D cube with the z-scale per pixel
        
        z3d = np.zeros(Tg.size, dtype='float32')
        for ii in range(z.size):
            z3d[ii,:,:] = z[ii]
            
    #
    # Now we can create simple Fits cube that the depth-optimizer can handle
    #
    nz, ny, nx = Tg.shape
    nLevel_H = 6
    NVAR_TOT = 10 + nLevel_H
    
    cub = np.zeros((ny, nx, NVAR_TOT, nz), dtype='float32')
    cub[:,:,0,:] = z
    #cub[:,:,1,:] = ltau # Not needed, will be filled by the code on output
    cub[:,:,2,:] = Tg
    cub[:,:,3,:] = Vz
    #cub[:,:,4,:] = Vturb # If you 
    cub[:,:,5,:] = Bx
    cub[:,:,6,:] = By
    cub[:,:,7,:] = Bz
    cub[:,:,8,:] = Ne
    #cub[:,:,9::,:] = nH.transpose((1,2,0,3)) # Only if we have H from the simulation, otherwise leave empty

    #
    # write to disk
    #
    writefits('tmp.fits', cub)
    cub = 0 # Free memory
    
    #
    # execute the c++ code
    #
    nthreads        = 63
    convert_units   = 0
    smooth_opt      = 15
    temperature_cut = 50000
    ltau_max        = 2
    new_ndep        = nz

    of.system('./dephOpt.x tmp.fits tmp_opt.fits {0} {1} {2} {3} {4} {5}'.format(nthreads, convert_units, smooth_opt, temperature_cut, ltau_max, new_ndep))


    
    #
    # Then readback the optimized atmosphere
    #

    odir = 'synthesis/'
    os.system('mkdir '+odir)
    
    print("reading...")
    cub = np.ascontiguousarray(readfits('tmp_opt.fits'))
    print("done ...")

    #
    # Convert to SI
    #
    m[:,:,0] *= 1/M_TO_CM
    m[:,:,3] *= 1/M_TO_CM
    m[:,:,4] *= 1/M_TO_CM
    m[:,:,5] *= 1/TESLA_TO_GAUSS
    m[:,:,6] *= 1/TESLA_TO_GAUSS
    m[:,:,7] *= 1/TESLA_TO_GAUSS
    # In 64 bits for densities
    Ne        = np.float64(m[:,:,8,:].squeeze())*M3_TO_CM3
    rho       = np.float64(m[:,:,9,:].squeeze())/RHO_TRANS
    NH        = np.float64(m[:,:,10::,:].squeeze().transpose((2,0,1,3)))*M3_TO_CM3
    
    rh.make_xarray_atmos(odir+'modelin_slanted.h5', m[:,:,2,:].squeeze(), \
                         m[:,:,3,:].squeeze(), m[:,:,0,:].squeeze(), rho=rho, \
                         nH=NH, ne=Ne, Bz = m[:,:,7].squeeze(), By=m[:,:,6].squeeze(), \
                         Bx=m[:,:,5].squeeze())

    # User defined wavelength array
    cw = 8542.091
    Dw = 1.75
    dw = 0.025
    ntw = int(2*Dw/dw +1)
    w = (np.arange(ntw)-ntw//2)*dw + cw

    rh.make_wave_file(odir+'wavelength_table.xr', new_wave=w*0.1, air=True)
