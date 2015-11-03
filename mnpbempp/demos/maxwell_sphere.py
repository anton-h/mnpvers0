# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 12:28:34 2015

@author: anton
"""
import sys
sys.path.append('/data/tho/my_bempp')
sys.path.append('/data/tho/bemppinst/lib/python2.7/site-packages')
# maxwell solution of dielectric sphere
import bempp.api
import bempp.api.operators.far_field
import numpy as np
import mnpbempp.bemsolver.bemret
import mnpbempp.particle
reload(mnpbempp.bemsolver.bemret)
# options for bempp
bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'

# energy of planewave excitation (eV)
ene = np.arange(1, 3.51, 0.1 )
# wavelength in nm
enei = 1240.0/ene
# number of wavelength values
n = enei.shape[0]

# make table of dielecric functions
import mnpbempp.material.material as material
epsd = material.epsdrude( 'silver' )
epsc = material.epsconst( 1 )
# store sourounding medium in first entry of list ( needed by some routines )
epstab = [ epsc, epsd ]

# make array for scattering data
ext = np.array(enei, float)

# import grid from file
grid = bempp.api.import_grid('/data/tho/disk2.msh')
# create comparticle object
p = mnpbempp.particle.comparticle( epstab, grid )

# create bemsolver for particle
bem = mnpbempp.bemsolver.bemret.bemret( p )

# create excitation   
polarisation = np.array([1,0,0])
direction = np.array([0,0,1])
import mnpbempp.simulation.planewave.excitation
exc = mnpbempp.simulation.planewave.excitation.exc( polarisation, direction) 

# loop over energies
for j in range(0,n):
  j  
  # solve with given excitation and wavelength
  sig = bem.solve( enei[j], exc( p, enei[j] ) )
  # compute extinction
  ext[j] = exc.ext( sig )
## plot solution and compare with mie solution  
from pymiecoated import Mie
#radius of sphere
R = 75
xs = 2*np.pi*R/enei

ext_mie = np.array(ene, float)
for imie in range(n):

    mie = Mie(x=xs[imie],eps=epstab[1](enei[imie]))
    ext_mie[imie]=mie.qext()
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#%matplotlib qt
#plt.plot( ene, ext,'ro', ene, ext_mie , 'b' )
plt.plot(ene,ext/(1239.8*8*np.pi),'ro',ene,ext_mie,'b')
