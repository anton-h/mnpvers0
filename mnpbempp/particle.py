# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 09:24:39 2015

@author: tho
"""
import bempp.api

class comparticle:

  def __init__( self, epstab, grid ):
    self.epstab = epstab
    self.grid = grid
    # create Raviart-Thomas space on grid
    self.space = bempp.api.function_space( grid, "RT", 0 )
   