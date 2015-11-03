# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:05:33 2015

@author: tho
"""  
# import needed modules for all classes
import numpy as np
from mnpbempp.misc import units
from mnpbempp import mnpbempp_dir

## dielectric constant
class epsconst:
  def __init__(self, dielectric_constant ):
    self.epsc = dielectric_constant
  def __call__(self, enei ):
    eps = self.epsc + 0*enei
    return eps
  
## drude dielectric function for different material parameters
class epsdrude:
# init  
  def __init__( self, material ):
    self.material = material
    #  atomic units
    hartree = 27.2116             #  2 * Rydberg in eV
    tunit = 0.66 / hartree        #  time unit in fs
    if  self.material  in ( 'Au', 'gold' ):
        self.rs = 3                     #  electron gas parameter
        self.eps0 = 10              #  background dielectric constant
        self.gammad = tunit / 10        #  Drude relaxation rate
    elif self.material in ( 'Ag', 'silver' ):
        self.rs = 3
        self.eps0 = 3.3
        self.gammad = tunit / 30
    elif self.material in ( 'Al', 'aluminum' ):
        self.rs = 2.07
        self.eps0 = 1
        self.gammad = 1.06 / hartree
    else:
        raise IOError( 'Material name unknown' )
    #  density in atomic units
    density = 3 / ( 4 * np.pi * self.rs ** 3 )
    #  plasmon energy
    self.wp = hartree * np.sqrt( 4 * np.pi * density )
    #  save values
    self.gammad = self.gammad * hartree
# methods    
  def __call__( self, enei ):
    #  convert to eV
    w = units.eV2nm / enei
    #  dielectric function and wavevector
    eps = self.eps0 - self.wp ** 2 / ( w * ( w + 1j * self.gammad ) );
#    #  wavenumber
#    k = 2 * np.pi / enei * np.sqrt( result )
    return eps

## tabulated dielectric numbers from experiment    
class epstable:
  # import needed modules for class table  
  def __init__( self, material ):
    from scipy.interpolate import interp1d
    # import tabulated data
    data = np.loadtxt(
          '{0}{1}{2}{3}'.format( mnpbempp_dir.path_mnpbempp,'/mnpbempp/material/',
           material, '.dat' ), skiprows = 2 )            
    self.enei_tab = units.eV2nm / data[ :, 0 ]
    self.n_tab    =               data[ :, 1 ]
    self.k_tab    =               data[ :, 2 ]
    self.n_interp = interp1d( self.enei_tab, self.n_tab )
    self.k_interp = interp1d( self.enei_tab, self.k_tab )
    
  def __call__( self, enei ):
    self.enei = enei
    #  assert that energy is in range
    assert( np.logical_and( np.min( self.enei_tab ) <= np.min(     self.enei ),
            np.max( self.enei     ) <= np.max( self.enei_tab ) ) ) 
    #  real and imaginary part of refractive index
    ni = self.n_interp( self.enei )
    ki = self.k_interp( self.enei )
    #  dielectric function
    eps = ( ni + 1j * ki ) ** 2
#    %  wavenumber
#    k = 2 * pi ./ enei .* sqrt( eps );
    return eps
  