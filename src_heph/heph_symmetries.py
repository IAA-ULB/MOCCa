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
# (i)  J. Dobaczewski et al. Phys. Rev. C  62, 014310 (2000). 
# (ii) J. Dobaczewski et al. Phys. Rev. C  62, 014311 (2000).
#
# The group D^T_2h is denoted as
# 
#  { 1,P,R_x,R_y,R_z,S_x,S_y,S_z,T,P^T,R_x^T,R_y^T,R_z^T,S_x^T,S_y^T,S_z^T}
#  (Eq. 3 in (ii))
# 
# and I will use this ordering to encode symmetry information in the .wf files.
#
#-------------------------------------------------------------------------------
import itertools
#-------------------------------------------------------------------------------
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
#  acting ONLY on the spatial coordinates, i.e. with per = [+1,+2,+3,+4] is 
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
  def __init__(self, perm, coord, linear, hermitian):
    # These are signed(!) permutations
    self.permutation = perm
    self.coord     = coord
    self.linear    = linear 
    self.hermitian = hermitian

  def __eq__(self, other):
    c = True

    if(self.permutation != other.permutation):
      c = False
    if(self.coord       != other.coord):
      c = False
    if(self.linear      != other.linear):
      c = False
    if(self.hermitian   != other.hermitian):
      c = False

    return c

  def __ne__(self, other):
    return not self.__eq__(other)


  def __hash__(self):
    # Explicitly make these things hashable
    return hash((self.permutation, self.coord, self.linear, self.hermitian))


  def __repr__(self):
    return "Symmetry"

  def __str__(self):
    return "S= [%+1d, %+1d, %+1d, %+1d]"%(self.permutation[0], \
                                          self.permutation[1], \
                                          self.permutation[2], \
                                          self.permutation[3]) \
        + ' x=%+1d,y=%+1d,z=%+1d'%(self.coord[0], self.coord[1], self.coord[2])\
        + ' linear    %5s '%str(self.linear) \
        + ' hermitian %5s '%str(self.hermitian)

  def __lt__(self, other):
     """
      Comparison operator for sorting lists of symmetries.
     """
     c1 = 0
     c2 = 0
     for i in range(3):
      if(self.coord[i] < 0):
        c1 = c1+1
      if(other.coord[i] < 0):
        c2 = c2+1

     if(c1 == c2):
       t1 = (abs(self.permutation[0]) ,abs(self.permutation[1]) , \
             abs(self.permutation[2]) ,abs(self.permutation[3]))
       t2 = (abs(other.permutation[0]),abs(other.permutation[1]), \
             abs(other.permutation[2]),abs(other.permutation[3]))
       return t1 < t2
     else:         
       return c1<c2 

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

    coord = [a*b for a,b in zip(S1.coord,S2.coord)]
    C = symmetry(perm,coord, lin, her )
    return C

class symmetry_option():
  """
      Practical container class to transmit all relevant symmetry combinations
      between different routines.
  """
  def __init__(self, generators, syms, combs, ReduceAxes, desc):
    self.generators = generators
    self.syms       = syms
    self.combs      = combs
    self.ReduceAxes = ReduceAxes 
    self.desc       = desc
      
    # Then, we check if there is a "timelike" symmetry conserved, i.e. an 
    # antilinear, antihermitian one
    self.timelike = False
    for g in generators:
      if( not g.linear and not g.hermitian):
        self.timelike = True

def symmetryencoding(so):
    """
       We encode all the information on symmetry choices into a simple string of
       ones and zeros. This is not very human-readable, but it is easy to 
       transfer in different coding languages.

       The first sixteen numbers indicate which generators of the 
       single-particle group have been chosen, in this ordering

       1 P [s] R_x R_y R_z [s] S_x S_y S_z [s] T P^T [s] 
                                        R_x^T R_y^T R_z^T [s] S_x^T S_y^T S_z^T

       where [s] indicates a space. Admittedly, the first entry in this list is 
       useless, but this way the group is at least complete.

       The next three numbers indicate the reduction of which axes have been 
       asked for. 
      
       So, EV8-mode for example would read
       
          0 1 001 000 10 000 010 111
          1 P R_i S_i TP RTi STi ax   

       To indicate 
          Parity       : first one
          z-signature  : second one
          time-reversal: third one
          y-timesimplex: fourth one
    
          calculation in 1/8th of the box: final 111 

    """

    encoded = '0 ' # start with the identity operator
    symdic  = populatesymmetries()
    direc   = ['x', 'y', 'z'] 

    # Parity
    n = 0
    for g in so.generators:
      if g == symdic['P']:
        n = 1
    encoded = encoded + '%d '%n
    
    # Signatures  
    for i in range(3):
      d = direc[i]
      n = 0
      for g in so.generators:
        if g == symdic['R%s'%d]:
          n = 1
      encoded = encoded + '%d'%n

    encoded = encoded + ' '
    # Simplexes      
    for i in range(3):
      d = direc[i]
      n = 0
      for g in so.generators:
        if g == symdic['S%s'%d]:
          n = 1
      encoded = encoded + '%d'%n
    encoded = encoded + ' '
    # Time-reversal
    n = 0
    for g in so.generators:
      if g == symdic['T']:
        n = 1
    encoded = encoded + '%d'%n

    # Time-parity
    n = 0
    for g in so.generators:
      if g == symdic['PT']:
        n = 1
    encoded = encoded + '%d'%n
    encoded = encoded + ' '

    # Time-signatures
    for i in range(3):
      d = direc[i]
      n = 0
      for g in so.generators:
        if g == symdic['RT%s'%d]:
          n = 1
      encoded = encoded + '%d'%n
    encoded = encoded + ' '

    # Time-simplexes
    for i in range(3):
      d = direc[i]
      n = 0
      for g in so.generators:
        if g == symdic['ST%s'%d]:
          n = 1
      encoded = encoded + '%d'%n
    encoded = encoded + ' '

    # Finally, reduction of axes
    for i in range(3):
      encoded = encoded + '%d'%so.ReduceAxes[i]   
    return encoded

def initsymmetries(SYMSTRING, REDUCE):
    """
      Parse the strings SYMSTRING and REDUCE and use them to make various 
      decisions based on the users symmetry choices. 
      
      These decisions are also printed to STDOUT. 
      - - - - - - - - - - - - - - - - - 

      This function returns:
      1. Generators : the list of generators of the symmetry group asked for
                      by the user.
      2. syms : the list of symmetry relations that can be used to reduce the 
                calculation to parts of the box.
      3. combs: the combinations of the generators that combine to the 
                symmetry relations. 
      4. ReduceAxes: a list of three numbers, indicating whether or not 
                     a given Cartesian axis is reduced.
    """
    symdic = populatesymmetries()

    # for printing purposes
    direc = ['x', 'y', 'z']   

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

    independent = not CheckIndependency(generators)
    com         = not CheckGenerators(generators)

    if(not com or not independent):
      print ("Trouble with your choice of symmetry generators.")
      exit()

    # Generating all possible symmetry relations
    symrel, combs = gen_symrelations(generators )
    
    #---------------------------------------------------------------------------
    # We identify the axis-reduction asked for, and generate practical 
    # symmetry relations for functions on the mesh.
    #---------------------------------------------------------------------------
    ReduceAxes = [0,0,0]
    for i in range(3):
      ReduceAxes[i] = int(REDUCE[i])

    c = 0
    for g in generators:
      if( g.linear or g.hermitian ):
        c = c + 1
    if(c > sum(ReduceAxes)):
      print (" The chosen symmetries allow for the reduction of more axes.")  
      print ("  Stopping.")
      exit()
    elif(c < sum(ReduceAxes)):
      print (" The chosen symmetries do not allow for this reduction.")
      print ("  Stopping.")
      exit()  
    syms, combs =  choose_symrelations(symrel, combs, ReduceAxes)
    # Some sanity checks
    for i in range(3):
      if(ReduceAxes[i] ==1 and (not syms[i])):
        print ("  Could not identify a relation to reduce axis %s."%direc[i])
        print ("  Stopping.")
        exit()

    so = symmetry_option(generators, syms, combs, ReduceAxes, SYMSTRING)
    return so

def printsymmetryoption(so):
    """
      Print the properties of this set of symmetry options.
    """

    # for printing purposes
    direc = ['x', 'y', 'z']   

    symdic = populatesymmetries()
    st = ''
    for g in so.generators:
      for s in symdic.keys():
        if (g == symdic[s]):
          st = st + s + ','
      
    if(len(so.generators) >=1):
      st = st [:-1]

    independent = not CheckIndependency(so.generators)
    com         = not CheckGenerators(so.generators)
    print ("     Generators: ", st )
    print ("      => Correct commutation relations?: " , com )
    print ("      => All independent?              : " , independent)

    print ("     Total number of symmetry relations: ", len(so.syms))    
    print ("     Time-like symmetry                : ", so.timelike)
    print ("  Reduction asked for:")
    print ("   (x,y,z) = (%d, %d, %d)"\
                       %(so.ReduceAxes[0], so.ReduceAxes[1], so.ReduceAxes[2]))

    for i in range(3):
      st = ''
      if(so.ReduceAxes[i] == 1):
        for l in range(len(so.combs[i])):
          for key in symdic.keys():
            if(so.combs[i][l] == symdic[key]):
              st = st + key + ','

        print ("   ", direc[i], " =>  %10s, "%st[:-1], so.syms[i])

    print (" Encoding string : ", symmetryencoding(so))

def populatesymmetries():
  """
    Returns a dictionary with all "named" symmetries in the symmetry groups
    of interest. This is primary useful to parse the users input and put 
    names on the symmetry relations deduced automatically.
  """

  symdic = {}
  #   Identity
  iden      = symmetry([+1,+2,+3,+4],[+1,+1,+1], True, True)
  symdic['1']  = iden
  # - Identity
  niden     = symmetry([-1,-2,-3,-4],[+1,+1,+1], True, True)
  symdic['-1'] = niden
  #   [P   psi](x,y,z,sigma) =          psi  [-x,-y,-z, sigma]
  P         = symmetry([+1,+2,+3,+4],[-1,-1,-1], True, True)
  symdic['P']  = P
  #   [T   psi](x,y,z,sigma) =    sigma psi^*[+x,+y,+z,-sigma]
  T           = symmetry([+3,-4,-1,+2],[+1,+1,+1], False, False)
  symdic['T'] = T
  #   [R_x psi](x,y,z,sigma) = -i       psi  [ x,-y,-z,-sigma]
  Rx           = symmetry([+4,-3,+2,-1],[+1,-1,-1], True, False)
  symdic['Rx'] = Rx 
  #   [R_y psi](x,y,z,sigma) =    sigma psi  [-x, y,-z,-sigma]
  Ry           = symmetry([-3,-4,+1,+2],[-1,+1,-1], True, False)
  symdic['Ry'] = Ry
  #   [R_z psi](x,y,z,sigma) = -i sigma psi  [-x,-y,+z, sigma]
  Rz           = symmetry([+2,-1,-4,+3],[-1,-1,+1], True, False)
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

  # Time-signatures
  RTx  = multiply(T, Rx)
  RTy  = multiply(T, Ry)
  RTz  = multiply(T, Rz)

  symdic['RTx'] = RTx
  symdic['RTy'] = RTy
  symdic['RTz'] = RTz

  # Time-parity
  PT   = multiply(P,T)
  symdic['PT'] = PT
  
  return symdic
#-------------------------------------------------------------------------------
def gen_symrelations(gen):
  """
     Starting from a set of single-particle symmetry operators, this routine 
     identifies all the symmetry relations that could be used to simplify the 
     numerical representation. 

     For an operator S that gives rise to either eigenstates or invariants, we 
     have that 
      
      [S psi](x,y,z,sigma) = p psi(x,y,z,sigma)

     where p is a number of modulus one. Depending on whether p is real or 
     imaginary there is a difference in the final symmetry relations imposed on
     the independent components of the spwfs. 

  """

  # First, we remove all antilinear, antihermitian symmetry operators. These
  # do not give rise to symmetry relations for individual spwfs.
  effgen = []
  for g in gen:
    if(  g.linear or g.hermitian ) :
       effgen.append(g)

  # We need to return the actual symmetry relations (sg) but we also return 
  # the combination (combs) of symmetry operators that make up the symmetry 
  # relations.
  sg    = []
  combs = []

  # We will construct by brute-force all symmetry relations that could be used
  # to reduce the problem
  added = True
  k     = 0
  iden      = symmetry([+1,+2,+3,+4],[+1,+1,+1], True, True)
  
  while(added):
    k = k +1 

    added = False
    newlist = list(itertools.combinations(effgen,k))  

    for comb in newlist:
      sy = iden
      for l in range(k):
          sy = multiply(sy,comb[l])
        
      found = False
      for l in sg:
       if (l == sy):
         found = True

      if(not found):
        added = True
        sg.append(sy)
        combs.append(comb)

  # Sort the final symmetry group in a practical way
  sg, combs = zip(*sorted(zip(sg, combs)))

  return sg, combs

def choose_symrelations(symrel, combs, redu):   
  """
  Using all generated symmetry relations, we select the ones that are best 
  suited to our needs.

  (i)  We need a permutation that is [+-1,+-2,+-3,+-4] for linear operators
                                     [+-2,+-1,+-4,+-3] for antilinear operators

  (ii) We prefer symmetry relations that imply the least amount of coordinate
       reflections.

  """ 
  syms = [None, None, None] 
  selcombs = [None, None, None]
  used = [0] * len(symrel)

  cand =  {}
  axes = range(3)
  for i in axes:
    cand[i] = []
    if(redu[i] == 1):
     # enumerate all symmetry relations that involve the right axis
     for k,s in enumerate(symrel):
      if(s.coord[i] == -1 and used[k] == 0):
        if(      s.hermitian and abs(s.permutation[0]) == 1):
          cand[i].append(k)
        elif(not s.hermitian and abs(s.permutation[0]) == 2):
          cand[i].append(k)
  
  # Then we sort the axis on the number of candidate symmetries
  # I.e. do first the axis set with the smallest number of candidates
  axes = sorted(axes, key= lambda x: len(cand[i]))
  for i in axes:
     for j in range(len(cand[i])):
        if( used[cand[i][j]] == 0):
          used[cand[i][j]] = 1
          syms[i]          = symrel[cand[i][j]]
          selcombs[i]      = combs[cand[i][j]] 
          break
  return syms, selcombs


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
        niden     = symmetry([-1,-2,-3,-4],[+1,+1,+1], True, True)
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
