#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module that determines the symmetries imposed by Hephaestos on 
# Tantalus/Prometheus.
#-------------------------------------------------------------------------------
# We will use the definitions and notation of
#
# J. Dobaczewski et al. Phys. Rev. C  62, 014310 (2000). 
# J. Dobaczewski et al. Phys. Rev. C  62, 014311 (2000).
#
#-------------------------------------------------------------------------------
import itertools

# Whether or not signature is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
Signature     = 0
# Whether or not parity is a quantum number of the spwfs. 
#   0 Broken
#   1 Conserved
Parity        = 0
# Whether or not time-simplex is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
#
TimeSimplex   = 0
# Whether or not time-signature is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
#
TimeSignature = 0
# Whether or not time-reversal is conserved. Can be 0 (broken) or 1 (conserved).
TimeReversal  = 1


# Whether you want a reduction of each of the Cartesian axes. 
# Note that not all choices are compatible with all symmetry choices. 
# Hephaestos will complain if it feels something is amiss.
# Reduce_x/y/z = 1 if reduced, 0 if fully stored. 
ReduceAxes= [1,1,1]


#-------------------------------------------------------------------------------
#
# In a single-particle coordinate space with separated proton and neutron 
# states, a single-particle wavefunction is a complex function of space and
# spin variables:
#                   psi(x,y,z,sigma)
#
# All possible single-particle symmetry operators operate on this in particlar
# ways
#
#   [P   psi](x,y,z,sigma) =          psi  [-x,-y,-z, sigma]
#   [T   psi](x,y,z,sigma) =    sigma psi^*[+x,+y,+z,-sigma]
#   [R_x psi](x,y,z,sigma) = -i       psi  [ x,-y,-z,-sigma]
#   [R_y psi](x,y,z,sigma) =    sigma psi  [-x, y,-z,-sigma]
#   [R_z psi](x,y,z,sigma) = -i sigma psi  [-x,-y,+z, sigma]
#
# The action of the simplex, time-simplex, time-signature and time-parity 
# operators can be deduced from these actions. (We could also have left out
# one of the signature operators.)
#
# Representing the spin-part of the spwfs as
#
#  Psi = ( psi_1 + i psi_2 )
#        ( psi_3 + i psi_4 )
#
# All of the symmetry operators can be defined in terms of what they do as a
# 
#  (1) a signed permutation of the indices (1,2,3,4) 
#  (2) signs of the coordinates, i.e. (-1,1,1) for a sign in the x-coordinate
#
# For example, 
#
#    P  => ( (+1,+2,+3,+4), -1,-1,-1)
#    Rz => ( (+2,-1,-4,+3), -1,-1,+1)
#  
# Combination of symmetry elements can then be achieved through
#   (1) combination of permutation of indices
#   (2) multiplication of all other signs. 
#
# So, for example, 
#
#    Sz => ( (+4,+3,+2,-1), +1, +1, -1 )
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Now, we also keep track of the properties of the symmetries. 
# Are they
#
#   *) linear    or antilinear? 
#   *) hermitian or antihermitian?
# 
# There must be some way to deduce that from the signed permutations + the 
# signs of the Cartesian coordinates, but I haven't found one. 
# (For hermeticity, there is the trivial observation that any operator 
#  acting ONLY on the spatial coordinates, i.e. with per = [1,2,3,4] is 
#  hermitian.)
#  
# Hence, a full classification becomes
#
# ( permutation, +-1, +-1, +-1,  linear?, hermitian? )
#
# where linear and hermitian are logicals. 
#-------------------------------------------------------------------------------
class symmetry():
  """
    Practical definition of a symmetry in terms of a signed permutation 
    of indices and signs on the x,y,z coordinates
  """
  def __init__(self, perm, x,y,z, linear, hermitian):
    # These are signed(!) permutations
    self.permutation = perm
    self.x         = x
    self.y         = y
    self.z         = z
    self.linear    = linear 
    self.hermitian = hermitian

  def __eq__(self, other):
    c = True

    if(self.permutation != other.permutation):
      c = False
    if(self.x != other.x):
      c = False
    if(self.y != other.y):
      c = False
    if(self.z != other.z):
      c = False
    if(self.linear    != other.linear):
      c = False
    if(self.hermitian != other.hermitian):
      c = False

    return c

  def __ne__(self, other):
    return not self.__eq__(other)


  def __hash__(self):
    # Explicitly make these things hashable
    return hash((self.permutation, self.x, self.y, self.z, self.linear, self.hermitian))


  def __repr__(self):
    return "Symmetry"

  def __str__(self):
    return "S= [%+1d, %+1d, %+1d, %+1d]"%(self.permutation[0], \
                                          self.permutation[1], \
                                          self.permutation[2], \
                                          self.permutation[3]) \
              + ' x=%+1d,y=%+1d,z=%+1d'%(self.x, self.y, self.z)\
              + ' linear    %5s '%str(self.linear) \
              + ' hermitian %5s '%str(self.hermitian)

  def __lt__(self, other):
     """
     """
     t1 = (abs(self.permutation[0]),abs(self.permutation[1]),abs(self.permutation[2]),abs(self.permutation[3]))
     t2 = (abs(other.permutation[0]),abs(other.permutation[1]),abs(other.permutation[2]),abs(other.permutation[3]))

     if(t1 == t2):  
       c1 = 0
       if(self.x < 0):
        c1 = c1+1
       if(self.y < 0):
        c1 = c1+1
       if(self.z < 0):
        c1 = c1+1
    
       c2 = 0
       if(other.x < 0):
        c2 = c2+1
       if(other.y < 0):
        c2 = c2+1
       if(other.z < 0):
        c2 = c2+1

       if(c1 == c2):
         return self.permutation[0] < other.permutation[0]
       else:         
         return c1<c2 

     else:
      return t1 < t2
#     c1 = 0
#     if(self.x < 0):
#      c1 = c1+1
#     if(self.y < 0):
#      c1 = c1+1
#     if(self.z < 0):
#      c1 = c1+1
#  
#     c2 = 0
#     if(other.x < 0):
#      c2 = c2+1
#     if(other.y < 0):
#      c2 = c2+1
#     if(other.z < 0):
#      c2 = c2+1

#     if(c1 == c2):
#       return abs(self.permutation[0]) < abs(other.permutation[0])
#     else:         
#       return c1<c2 

def multiply(S1, S2):
    """
    Combine two symmetry operators into a third symmetry:
          C  = S1 * S2
    """
    perm  = [0,0,0,0]
    for i in range(4):
        ind     = int(abs(S1.permutation[i])-1)
        s       = S1.permutation[i]/abs(S1.permutation[i])
        perm[i] = s*S2.permutation[ind]

    lin = S1.linear    == S2.linear
    her = S1.hermitian == S2.hermitian
    C = symmetry(perm,S1.x*S2.x,S1.y*S2.y,S1.z*S2.z, lin, her )
    return C

#-------------------------------------------------------------------------------
# Populate all the symmetries with correct behavior
#-------------------------------------------------------------------------------
symdic = {}
#   Identity
iden      = symmetry([+1,+2,+3,+4],+1,+1,+1, True, True)
symdic['1']  = iden
# - Identity
niden     = symmetry([-1,-2,-3,-4],+1,+1,+1, True, True)
symdic['-1'] = niden
# Multiplication with i
imag      = symmetry([-2,+1,-4,+3],+1,+1,+1, True, False)
symdic['i'] = imag
#   [P   psi](x,y,z,sigma) =          psi  [-x,-y,-z, sigma]
P         = symmetry([+1,+2,+3,+4],-1,-1,-1, True, True)
symdic['P']  = P

#   [T   psi](x,y,z,sigma) =    sigma psi^*[+x,+y,+z,-sigma]
T           = symmetry([-3,+4,+1,-2],+1,+1,+1, False, False)
symdic['T'] = T

#   [R_x psi](x,y,z,sigma) = -i       psi  [ x,-y,-z,-sigma]
Rx           = symmetry([+4,-3,+2,-1],+1,-1,-1, True, False)
symdic['Rx'] = Rx 

#   [R_y psi](x,y,z,sigma) =    sigma psi  [-x, y,-z,-sigma]
Ry           = symmetry([-3,-4,+1,+2],-1,+1,-1, True, False)
symdic['Ry'] = Ry

#   [R_z psi](x,y,z,sigma) = -i sigma psi  [-x,-y,+z, sigma]
Rz           = symmetry([+2,-1,-4,+3],-1,-1,+1, True, False)
symdic['Rz'] = Rz

# Simplexes
Sx   = multiply(P, Rx)
Sy   = multiply(P, Ry)
Sz   = multiply(P, Rz)

symdic['Sx'] = Sx
symdic['Sy'] = Sy
symdic['Sz'] = Sz

# Time-simplexes
STx  = multiply(T, Sx)
STy  = multiply(T, Sy)
STz  = multiply(T, Sz)

symdic['STx'] = STx
symdic['STy'] = STy
symdic['STz'] = STz

print ('STx', STx)

# Time-parity
PT   = multiply(P,T)
symdic['PT'] = PT


def populate_symgroup(gen):
  """
    Take the symmetry generators and form the whole group. 
    The chief aim of this routine is to obtain a set of independent symmetry    
    relations that can be employed to simplify the calculation.

    1. Generate all possible unique multiplications of the generators
    2. Then, convert these symmetry relations to a specific convention
        a. if( permutation[0]> permutation[1] )
           multiply with  i, such that permutation[0] < permutation[1]

    Note that 

    a) step 1. is accomplished using a simple brute-force enumeration method.
    b) the additional multiplication is usually needed to create relations 
       with the four components of the wavefunctions that are one-to-one, i.e.
        psi_1 is related to the symmetry-transformed psi_1.

            R_z psi = \pm i psi
       => i R_z psi =-\pm   psi 
       => psi_{1-4} =-\pm   psi


  """
  sg    = []
  combs = []
  added = True
  k     = 0
  while(added):
    k = k +1 

    added = False
    newlist = list(itertools.combinations(range(len(gen)),k))  

    for comb in newlist:
      sy = iden
      for l in range(k):
          sy = multiply(sy,gen[comb[l]])
        
      found = False
      for l in sg:
       if (l == sy):
         found = True

      if(not found):
        added = True
        sg.append(sy)
        combs.append(comb)

  for i in range(len(sg)):
    if(abs(sg[i].permutation[0])>abs(sg[i].permutation[1])):
      sg[i] = multiply(imag,sg[i])
      t = list(combs[i])
      t.append(len(gen))
      combs[i] = tuple(t)
#  for i in range(len(sg)):
#    if( sg[i].permutation[0]<0):
#      sg[i] = multiply(niden,sg[i])
#      t = list(combs[i])
#      t.append(len(gen)+1)
#      combs[i] = tuple(t)

  # Sort the final symmetry group in a practical way
  sg, combs = zip(*sorted(zip(sg, combs)))
  return sg, combs


#-------------------------------------------------------------------------------

def initsymmetries(SYMSTRING, REDUCE):
    #---------------------------------------------------------------------------
    # Parse the SYMSTRING and use it to initialize the whole module.      
    #---------------------------------------------------------------------------
    # Parse the string 
    global Signature, Simplex, Parity,TimeSimplex,TimeReversal,TimeSignature
    global ReduceAxes    

    syms       = SYMSTRING.split(',')
    generators = []
    for s in syms:
      found = False
      for key in symdic.keys():
        if s == key:
            generators.append(symdic[key])
            found = True

      if(not found):
        print('Error parsing symmetries. %s is not recognized.'%s)


    if(len(generators)>4):
      print ('Maximum length of generators of the symmetry group = 4.')
      exit()

    print ("  Symmetry information" )
    st = ''
    for s in syms:
      st = st + s + ','
    if(len(syms) >=1):
      st = st [:-1]

    independent = not CheckIndependency(generators)
    com         = not CheckGenerators(generators)
    print ("     Generators: ", st )
    print ("      => Correct commutation relations?: " , com )
    print ("      => All independent?              : " , independent)

    if(not com or not independent):
      print ("Trouble with your choice of symmetry generators.")
      exit()
    # Generating the full symmetry group  
    generators = generators
    symgroup, combs = populate_symgroup(generators )
    
    print ("     Total size of the symmetry group  : ", len(symgroup))    

    for i,s in enumerate(symgroup):
##      if(abs(s.permutation[0]) == 1):
        print (i,s)
    #---------------------------------------------------------------------------
    # We identify the axis-reduction asked for, and generate practical 
    # symmetry relations for functions on the mesh.
    #---------------------------------------------------------------------------
    print ("  Reduction asked for:")
      
    redu_x = int(REDUCE[0])
    redu_y = int(REDUCE[1])
    redu_z = int(REDUCE[2])
    print ("   (x,y,z) = (%d, %d, %d)"%(redu_x, redu_y, redu_z))

    c = 0
    for g in generators: # Do not count the imaginary unit at the end
      if( g.linear or g.hermitian ):
        c = c + 1
    if(c > redu_x + redu_y + redu_z):
      print (" The chosen symmetries allow for the reduction of more axes.")  
      print ("  Stopping.")
      exit()
    elif(c < redu_x + redu_y + redu_z):
      print (" The chosen symmetries do not allow for this reduction.")
      print ("  Stopping.")
      exit()  


    # Finding the practical symmetry relations that will allow us to reduce the 
    # axes. 
    xsym, xcomb,ysym, ycomb,zsym,zcomb = GenSymRelations(symgroup, combs, redu_x, redu_y, redu_z)

    # Some sanity checks
    if(redu_x ==1 and (not xsym)):
      print ("  Could not identify a symmetry relation to reduce the X-axis.")
      print ("  Stopping.")
      exit()
    if(redu_y ==1 and (not ysym)):
      print ("  Could not identify a symmetry relation to reduce the Y-axis.")
      print ("  Stopping.")
      exit()
    if(redu_z ==1 and (not zsym)):
      print ("  Could not identify a symmetry relation to reduce the Z-axis.")
      print ("  Stopping.")
      exit()
    
    generators = generators + [imag, niden]
    if(redu_x == 1):
      st = ''
      for i in xcomb:
        for a in symdic.keys():
          if(generators[i] == symdic[a]):
            st = st + a + ','

      print ("  ", xsym, st)
    if(redu_y == 1):
      st = ''
      for i in ycomb:
        for a in symdic.keys():
          if(generators[i] == symdic[a]):
            st = st + a + ','
      print ("  ", ysym, st)
    if(redu_z == 1):
      st = ''
      for i in zcomb:
        for a in symdic.keys():
          if(generators[i] == symdic[a]):
            st = st + a + ','
      print ("  ", zsym, st)


def CheckGenerators(gen):
  """
    Small routine to check if the generators commute/anti-commute as they 
    should to obtain a set of symmetry-operators that can impose the right 
    structure on the spwfs.
  """

  problem = False
  for f in gen:
    for g in gen:
      fg = multiply(f,g)
      gf = multiply(g,f)

      if(f.linear and g.linear):
        # Both linear operators, textbook case
        if(fg != gf):
           problem =  True

      elif( f.hermitian and g.hermitian):
        # Both hermitian; one is linear and one is antilinear
        if(fg != gf):
           problem =  True
      elif( (not f.hermitian) and (not f.linear)):
        if(fg != gf):
           problem =  True
      elif( (not g.hermitian) and (not g.linear)):
        if(fg != gf):
           problem =  True
      else:
        # One is antihermitian, linear, the other antihermitian, linear.
        # They need to ANTICOMMUTE in this case
        gf = multiply(niden, gf)
        if(fg != gf):
           problem =  True
  return problem     

def CheckIndependency(gen):
  """
    Check whether or not the set of generators are independent, i.e. there is
    not one that can be generated by multiplying the others. 
  """

  problem = False
  for i in range(len(gen)): 
    nl = gen[i+1:] + gen[0:i]
      
    for c in nl:
      if (c == gen[i]):
        problem  = True

    # all two-combinations
    two = list(itertools.combinations(nl,2))  
    for cb in two:
      c = multiply(cb[0], cb[1])
      if (c == gen[i]):
        problem  = True

    # all three-combinations
    three = list(itertools.combinations(nl,3))  
    for cb in three:
      c = multiply(cb[0], cb[1])
      c = multiply(c, cb[2])
      if (c == gen[i]):
        problem  = True

  return problem

def GenSymRelations(symgroup, combs, redu_x, redu_y, redu_z):   
  """
  Using the full symmetry-group, we try to find practical symmetry relations  
  to use to reduce the calculation to smaller meshes.

  Note that there is no real "preference" among relations that is built into
  this routine. However, if the input symgroup is sorted (according to the 
  __lt__ operator defined above), then this routine will return symmetry 
  relations that employ the least of "reflections" possible. I.e. it will
  prefer symmetry relations with X coordinate signs over other ones with 
  Y coordinate changes if X < Y.
  """ 

  xsym = None 
  xcomb = None
  ysym = None 
  ycomb = None
  zsym = None
  zcomb = None

  used = [0] * len(symgroup)

  if(redu_x == 1):
   for i,s in enumerate(symgroup):
    if(s.x == -1):
        xsym = s
        xcomb= combs[i]
        used[i] = 1
        break      

  if(redu_y == 1):
   for i,s in enumerate(symgroup):  
      if(used[i] == 1):
        continue
      if(s.y == -1):
          ysym = s   
          ycomb= combs[i]
          used[i] = 1
          break      
  
  if(redu_z == 1):
   for i,s in enumerate(symgroup):  
      if(used[i] == 1):
        continue
      if(s.z == -1):
          zsym = s
          zcomb= combs[i]
          used[i] = 1
          break      
  return xsym, xcomb, ysym, ycomb, zsym, zcomb

#def CheckSymmetryChoice() :
#    #---------------------------------------------------------------------------
#    # Checking whether the choices made for the symmetries are valid.
#    # The options checked for are trivial for now.
#    #---------------------------------------------------------------------------
#    global Signature, Simplex, Parity, TimeSimplex, TimeReversal, TimeSignature
#    
#    # First off, does anything take illegal values 
#    if(Signature > 3 or Signature < 0 ):
#            print ('Illegal value for signature.')
#            exit()
#    if(Simplex > 3 or Simplex < 0 ):
#            print ('Illegal value for simplex.')
#            exit()
#    if(TimeSimplex > 3 or TimeSimplex < 0 ):
#            print ('Illegal value for time-simplex.')
#            exit()
#    if(TimeSignature > 3 or TimeSignature < 0 ):
#            print ('Illegal value for time-signature.')
#            exit()
#    if(Parity > 1 or Parity < 0 ):
#            print ('Illegal value for parity.')
#            exit()
#    if(TimeReversal > 1 or TimeReversal < 0 ):
#            print ('Illegal value for time-reversal.')
#            exit()

#def CheckReduction():
#    #---------------------------------------------------------------------------
#    # Check whether the asked-for reduction of axes is valid, given the 
#    # particular combination of symmetries asked for.
#    #---------------------------------------------------------------------------

#    global Signature, Simplex, Parity, TimeSimplex, TimeReversal, TimeSignature
#    global ReduceAxes
#    
#    # Check whether we are not over-or under-using the symmetry information.
#    s = sum(ReduceAxes)
#    p = 0
#    if(Signature != 0):
#             p = p + 1 
#    if(Simplex != 0):
#             p = p + 1 
#    if(TimeSimplex != 0):
#             p = p + 1 
#    if(TimeSignature != 0):
#             p = p + 1 
#    if(Parity != 0):
#             p = p + 1 

#    if(s > p) :
#     print ('You asked for more reduction of axes than the symmetries can cover.')
#     exit()
#    elif(s < p):
#     print ('The symmetries you demanded provide more reduction of axes ' \
#         + 'than you asked for.')
#     exit()

#    #---------------------------------------------------------------------------
#    # Check more in detail
#    for i in range(3):
#      Valid = False      
#      if(ReduceAxes[i] == 1):
#        if(Simplex == i +1):
#            Valid = True
#        if(TimeSimplex == i +1):        
#            Valid = True
#        if(Signature != i +1 and Signature != 0):        
#            Valid = True
#        if(Parity != i +1 and Parity != 0):        
#            Valid = True
#        if( not Valid) :
#            print 'The symmetry choice does not support reducing axis %d'%i
#            exit

#    #---------------------------------------------------------------------------


#def PrintSymmetry():
#        # Outputs all of the information on symmetry conservation.
#        global Signature,Simplex,Parity,TimeSimplex,TimeReversal,TimeSignature
#        global ReduceAxes 
#        
#        words      = ['Broken', 'Conserved']
#        directions = ['x', 'y', 'z']
#        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
#        print ' Explicit spwf symmetries '
#        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
#             
#        print ' Time-reversal: ' + words[TimeReversal]
#        print ' Parity       : ' + words[Parity]

#        print '                 X Y Z '
#        print ' Signature     : %d %d %d'%(Signature    == 1, Signature     == 2, Signature     == 3)
#        print ' Simplex       : %d %d %d'%(Simplex      == 1, Simplex       == 2, Simplex       == 3)
#        print ' TimeSignature : %d %d %d'%(TimeSignature== 1, TimeSignature == 2, TimeSignature == 3)                
#        print ' TimeSimplex   : %d %d %d'%(TimeSimplex  == 1, TimeSimplex   == 2, TimeSimplex   == 3)
#        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
#        print         
#        for i in range(3):
#                if ReduceAxes[i] == 1 :
#                        print ' %s-axis : positive half stored'%directions[i]
#                else:
#                        print ' %s-axis : fully         stored'%directions[i]
#        print              
#        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

