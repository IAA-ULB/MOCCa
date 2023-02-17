import numpy         as np
from read_density    import *
from scipy.special   import sph_harm
from scipy.integrate import romb
import matplotlib.pyplot  as plt

def decompose_density_axial(fname, rr, ls=[0], power=5):
  """
    Obtain the multipole components of an (axially symmetric) density, i.e.    
    get some of the coefficients in 
    
    rho(r,theta) = sum_{lm} r_lm(r)  Re[Y_lm(theta,phi)]

    by explicitly integrating 
      
    r_lm(r) = int dtheta dphi

  """

  #--------------------------------------------
  # Read the density from fname
  nx, ny, nz, dx, den = read_den_file(fname)
  #--------------------------------------------
  # define a mesh in r
  den_r = np.zeros((len(rr), len(ls)))

  # Loop over multipole moments
  for k in range(len(ls)):
    l = ls[k]
    
    for i,r in enumerate(rr):

       #------------------------------------------------------------------------       
       # Now, we define a function of (theta, phi) that we want to integrate over
       def integrand(theta, phi=0):
          
          # Construct the 3D coordinate of the point at (r,theta,phi)
          try:
            point = np.zeros((len(theta),3))
          except TypeError:
            point = np.zeros((1,3))
          point[:,0] = r*np.sin(theta)*np.cos(phi)  # x = r sin(theta) cos(phi)
          point[:,1] = r*np.sin(theta)*np.sin(phi)  # y = r sin(theta) sin(phi)
          point[:,2] = r*np.cos(theta)              # Z = r cos(theta)
       
          # Interpolate with Lagrange functions
          inter_den = Interpolate_function(point,den[:,0:3],den[:,4],+1,+1,+1)
       
          # multiply with the (real part of the) spherical harmonic
          Ylm = np.real(sph_harm(0,l,0,theta)) 
          return inter_den *Ylm * np.sin(theta) 
          # sin(theta) is part of the surface element of integration

       # We don't take 0 itself, blows up for l > 0
       # Integration from 0->pi/2
       theta_r  = np.linspace(0,np.pi/2,2**power+1,endpoint=True)     
       dtheta = theta_r[1]-theta_r[0]
       
       integ = integrand(theta_r)
       den_r[i,k] = romb(integ, dx=dtheta)              # Romberg integration
       den_r[i,k] = den_r[i,k]              * 2 * np.pi # Integration over phi
       den_r[i,k] = den_r[i,k]              * 2         # Full theta integration
    
  return den_r


