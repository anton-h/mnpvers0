# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:28:48 2015

@author: tho
"""
## only x polarisation an wavevec in z-direction ...
import bempp.api
import bempp.api.operators.far_field
import numpy as np

class exc:
  def __init__( self, polarisation, direction ):
    ## only x polarisation and wavevec in z-direction ...
    if  ( any( polarisation != np.array([1,0, 0 ]) ) | any( direction != np.array([0, 0, 1]) ) ):
      raise IOError( 'Till now, only one polarisation and direction allowed' )

    self.pol = polarisation
    self.dir = direction
  
  # evaluate Dirichlet and Neumann trace at grid 
  def __call__( self, p, enei ):
    self.p = p 
    self.kext = 2*np.pi / enei * np.sqrt( p.epstab[0](enei) )
      # incoming field  
    def evalIncField(point):
      x, y, z = point
      return np.array( [ np.exp( 1j * self.kext * z ), y * 0., z * 0. ] )
    
    ## Boundary conditions  return field    
    def evalincDT(point, normal, domain_index, result):
      field = evalIncField(point)
      result[:] = np.cross( field, normal, axis=0 )
#    return result

    def evalincNT( point, normal, domain_index, result ):
      x, y, z = point
      curl = np.array([x * 0., 1j * self.kext * np.exp(1j * self.kext * z), x * 0.])
      result[:] = np.cross(curl / (1j * self.kext), normal, axis=0)

    # create maxwell identity operator
    ido = bempp.api.operators.boundary.sparse.maxwell_identity( p.space )
  
    self.DTinc = bempp.api.GridFunction( p.space, fun=evalincDT )
    self.NTinc = bempp.api.GridFunction( p.space, fun=evalincNT )

    # right hand side
    rhs = np.hstack( ( ( ido*self.NTinc ).projections(), ( ido*self.DTinc ).projections() ) )
    return rhs
  
  # compute extinction for given polarisation and direction  
  def ext( self, sig ): 
    DText = sig[0] - self.DTinc
    NText = sig[1] - self.NTinc

    # Create the necessary potential operators
    
    slfPot = bempp.api.operators.far_field.maxwell.electric_field( self.p.space, np.vstack( self.dir ), self.kext)
    dlfPot = bempp.api.operators.far_field.maxwell.magnetic_field( self.p.space, np.vstack( self.dir ), self.kext)
    ffp    = slfPot.evaluate(NText) + dlfPot.evaluate(DText)
    Cext   = 4*np.pi / (self.kext ** 2 ) * np.real( 
                              1j * self.kext * np.inner( ffp[:,0], self.pol ) )
    return Cext
  