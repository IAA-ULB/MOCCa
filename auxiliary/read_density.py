import numpy as np
import matplotlib.pyplot as plt

def Lagrange_interpolation_function(x, xr, nx, dx) :
   
  """ 
    Calculate the value of a (one-D) Lagrange interpolation function 
    (associated with a reference point x_r) at a (different) point x:
                 
      f_r(x) = (2 nx)^{-1} sin[ a * (x - x_r) ] / [  sin (a * (x - x_r)/(2nx)]

      [Eq. (16)  in W. Ryssens et al., PRC 064318 (2015) ]
    
    where 
        a   = pi/dx
        x_r = coordinate of the reference point
        dx  = spacing of the mesh
        nx  = number of points in the given direction 
              (in a symmetry-unrestricted box)
  """  

  a = np.pi/dx
  N = 1.0/(2*nx)

  f      = N * np.sin(a*( x - xr))/(np.sin(a *N* (x - xr)))

  # We have to take care that (x - xr) is not too close to zero, the code
  # will throw us a NaN instead of the correct answer (1)
  ind = np.where(abs(x - xr) < 1e-10 )    
  f[ind] = 1.0

  return f

def Interpolate_function(points,mesh_points,f_values,sx,sy,sz):
  """
    Interpolate a function given on a Cartesian Lagrange mesh to a new 
    set of points (which do not have to be from another Lagrange mesh).

    Input : 
      points     : numpy array containing new points
      mesh_points: numpy array containing the values of the Lagrange mesh points

      sx, sy, sz : signs of the transformation of the function under spatial 
                   reflections in every direction.
                   Should be +1 or -1 if there is a symmetry in said direction, 
                   0 if the symmetry is broken
  
    Output:
      values :  interpolated function values on the new points.

    Intermediate things:
      nx, ny, nz : number of mesh points in each direction.
      dx         : mesh spacing
  """

  dx = mesh_points[1,0] - mesh_points[0,0]

  assert (abs(sx)<=1)  
  assert (abs(sy)<=1)  
  assert (abs(sz)<=1)  
  
  # Note that in both cases (symmetry conserved or not) this returns 
  # HALF the number of mesh points in a full box.   
  nx = int(np.max(mesh_points[:,0])/dx + 0.5)
  ny = int(np.max(mesh_points[:,1])/dx + 0.5)
  nz = int(np.max(mesh_points[:,2])/dx + 0.5)

  den = np.zeros(len(points))
  for i in range(len(mesh_points[:,0])):

    xi = mesh_points[i,0]
    yi = mesh_points[i,1]
    zi = mesh_points[i,2]

    Lx =           Lagrange_interpolation_function(points[:,0],  xi, 2*nx, dx) 
    if(sx != 0):
      Lx = Lx + sx*Lagrange_interpolation_function(points[:,0], -xi, 2*nx, dx)

    Ly =           Lagrange_interpolation_function(points[:,1],  yi, 2*ny, dx) 
    if(sy != 0):
      Ly = Ly + sy*Lagrange_interpolation_function(points[:,1], -yi, 2*ny, dx)

    Lz =           Lagrange_interpolation_function(points[:,2],  zi, 2*nz, dx)
    if(sz != 0):  
      Lz = Lz + sz*Lagrange_interpolation_function(points[:,2], -zi, 2*nz, dx)

    den = den + Lx * Ly * Lz * f_values[i]
  
  return  den


def read_den_file(fname):
  """
   Read the numerical data with numpy and deduce some extra info on the mesh.
  """
  dat = np.loadtxt(fname)

  dx = dat[1,0] - dat[0,0]
  nx = int(np.max(dat[:,0])/dx + 0.5)
  ny = int(np.max(dat[:,1])/dx + 0.5)
  nz = int(np.max(dat[:,2])/dx + 0.5)
  
  if(np.min(dat[:,2])< 0.0):
    nz = 2*nz    

  return nx,ny,nz,dx, dat
