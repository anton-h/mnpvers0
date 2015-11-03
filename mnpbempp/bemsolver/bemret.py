# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:11:16 2015

@author: tho
"""
import bempp.api
import numpy as np
from numpy.linalg import solve as sol
class bemret:
  def __init__( self, p ):
    
    self.space = p.space
    self.epstab = p.epstab

#  def __call__( self ):
  def calderon( self, enei ):
    self.kext = 2*np.pi / enei * np.sqrt( self.epstab[0](enei) )
    self.kint = 2*np.pi / enei * np.sqrt( self.epstab[1](enei) )
    rho = self.kint / self.kext
    # create maxwell single and double layer operators
    slpext = bempp.api.operators.boundary.maxwell.electric_field( 
                                self.space, self.space, self.space, self.kext )
    dlpext = bempp.api.operators.boundary.maxwell.magnetic_field(
                                self.space, self.space, self.space, self.kext )    
    slpint = bempp.api.operators.boundary.maxwell.electric_field(
                                self.space, self.space, self.space, self.kint )       
    dlpint = bempp.api.operators.boundary.maxwell.magnetic_field( 
                                self.space, self.space, self.space, self.kint )             
                                
    lhsOp = bempp.api.BlockedOperator(2,2)
    lhsOp[0,0] = -(slpext + rho * slpint)
    lhsOp[0,1] = dlpext + dlpint
    lhsOp[1,0] = dlpext + dlpint
    lhsOp[1,1] = slpext + (1. / rho) * slpint 
    # convert into matrix
    lhsOp = np.vstack( [ np.hstack( [lhsOp[0,0].weak_form().A, lhsOp[0,1].weak_form().A ] ),
                         np.hstack( [lhsOp[1,0].weak_form().A, lhsOp[1,1].weak_form().A ] ) ] )
    return lhsOp
    
  def solve( self, enei, rhs ):
    self.mat = self.calderon( enei )
    x = sol( self.mat, rhs )
    # extract solution into Dirichlet and Neumann trace
    DTproj = x[ 0 : x.shape[0]/2 ]
    NTproj = x[ x.shape[0]/2: x.shape[0] ]
    # create grid functions from projections
    DTsol = bempp.api.GridFunction( self.space, projections=DTproj) 
    NTsol = bempp.api.GridFunction( self.space, projections=NTproj)
    # merge to list and return
    return [ DTsol, NTsol ]
    