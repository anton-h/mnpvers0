# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 08:41:07 2015

@author: tho
"""
# add path
# test dielectric function
# import modules
import numpy as np
import matplotlib.pyplot as plt
from mnpbempp.material import material as material
reload(material)

# name of material
mat = 'silver'
# construct dielectric functions
epsc = material.epsconst( 0.85+0.9j )
epsd = material.epsdrude( mat )
epst = material.epstable( mat )

enei = np.linspace( 200, 1000, num = 100 )

# plot
plt.plot( enei, np.real( epsc( enei ) ),'.', enei, np.imag( epsc( enei ) ),'.' )
plt.figure()
plt.plot( enei, np.real( epsd( enei ) ), '.b', enei, np.real( epst( enei ) ), '.r' )
plt.figure()  
plt.plot( enei, np.imag( epsd( enei ) ), '.b', enei, np.imag( epst( enei ) ), '.r' )
