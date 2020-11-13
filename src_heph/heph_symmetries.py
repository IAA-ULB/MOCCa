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
#   [R_y psi](x,y,z,sigma) =    sigma psi  [-x, y,-z, sigma]
#   [R_z psi](x,y,z,sigma) = -i sigma psi  [-x,-y,+z,-sigma]
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
#    Rz => ( (+4,+3,+2,-1), -1,-1,+1)
#  
# Combination of symmetry elements can then be achieved through
#   (1) combination of permutation of indices
#   (2) multiplication of all other signs. 
#
# So, for example, 
#
#    Sz => ( (+4,+3,+2,-1), +1, +1, -1 )
#
#-------------------------------------------------------------------------------
#
#
#  Step 1: Parse the generators of the subgroup from the input
#  Step 2: Generate all possible symmetries in the group that are actually
#          generated. 
#
#  Questions: 
#    *) Can I see from this type of excercise what are the antilinear symmetries?
#    *) 
#
#-------------------------------------------------------------------------------

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

def initsymmetries(SYMSTRING, REDUCE):
    #---------------------------------------------------------------------------
    # Parse the SYMSTRING and use it to initialize the whole module      
    #---------------------------------------------------------------------------
    # Parse the string 
    global Signature, Simplex, Parity,TimeSimplex,TimeReversal,TimeSignature
    global ReduceAxes    

    direc = { "x" : 1, "y" : 2 , "z" : 3}
  
    Parity       = 0 ; Signature     = 0 ; TimeSimplex = 0
    TimeReversal = 0 ; TimeSignature = 0 ; Simplex     = 0

    syms = SYMSTRING.split(',')
    for s in syms:
      if s[0] == 'P': 
        Parity = 1

      elif s[0] == 'R':
        if(s[1] == 'T'):
          TimeSignature = direc[s[2]]
        else:
          Signature = direc[s[1]]

      elif s[0] == 'S':
        if(s[1] == 'T'):
          TimeSimplex = direc[s[2]]                
        else:
          Simplex = direc[s[1]]

      elif s[0] == 'T':
        TimeReversal = 1
    #---------------------------------------------------------------------------
    # Double check the options
    CheckSymmetryChoice()
    CheckReduction()

    #---------------------------------------------------------------------------
    # Print what was chosen
    PrintSymmetry()

def CheckSymmetryChoice() :
    #---------------------------------------------------------------------------
    # Checking whether the choices made for the symmetries are valid.
    # The options checked for are trivial for now.
    #---------------------------------------------------------------------------
    global Signature, Simplex, Parity, TimeSimplex, TimeReversal, TimeSignature
    
    # First off, does anything take illegal values 
    if(Signature > 3 or Signature < 0 ):
            print ('Illegal value for signature.')
            exit()
    if(Simplex > 3 or Simplex < 0 ):
            print ('Illegal value for simplex.')
            exit()
    if(TimeSimplex > 3 or TimeSimplex < 0 ):
            print ('Illegal value for time-simplex.')
            exit()
    if(TimeSignature > 3 or TimeSignature < 0 ):
            print ('Illegal value for time-signature.')
            exit()
    if(Parity > 1 or Parity < 0 ):
            print ('Illegal value for parity.')
            exit()
    if(TimeReversal > 1 or TimeReversal < 0 ):
            print ('Illegal value for time-reversal.')
            exit()

def CheckReduction():
    #---------------------------------------------------------------------------
    # Check whether the asked-for reduction of axes is valid, given the 
    # particular combination of symmetries asked for.
    #---------------------------------------------------------------------------

    global Signature, Simplex, Parity, TimeSimplex, TimeReversal, TimeSignature
    global ReduceAxes
    
    # Check whether we are not over-or under-using the symmetry information.
    s = sum(ReduceAxes)
    p = 0
    if(Signature != 0):
             p = p + 1 
    if(Simplex != 0):
             p = p + 1 
    if(TimeSimplex != 0):
             p = p + 1 
    if(TimeSignature != 0):
             p = p + 1 
    if(Parity != 0):
             p = p + 1 

    if(s > p) :
     print ('You asked for more reduction of axes than the symmetries can cover.')
     exit()
    elif(s < p):
     print ('The symmetries you demanded provide more reduction of axes ' \
         + 'than you asked for.')
     exit()

    #---------------------------------------------------------------------------
    # Check more in detail
    for i in range(3):
      Valid = False      
      if(ReduceAxes[i] == 1):
        if(Simplex == i +1):
            Valid = True
        if(TimeSimplex == i +1):        
            Valid = True
        if(Signature != i +1 and Signature != 0):        
            Valid = True
        if(Parity != i +1 and Parity != 0):        
            Valid = True
        if( not Valid) :
            print 'The symmetry choice does not support reducing axis %d'%i
            exit

    #---------------------------------------------------------------------------


def PrintSymmetry():
        # Outputs all of the information on symmetry conservation.
        global Signature,Simplex,Parity,TimeSimplex,TimeReversal,TimeSignature
        global ReduceAxes 
        
        words      = ['Broken', 'Conserved']
        directions = ['x', 'y', 'z']
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
        print ' Explicit spwf symmetries '
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
             
        print ' Time-reversal: ' + words[TimeReversal]
        print ' Parity       : ' + words[Parity]

        print '                 X Y Z '
        print ' Signature     : %d %d %d'%(Signature    == 1, Signature     == 2, Signature     == 3)
        print ' Simplex       : %d %d %d'%(Simplex      == 1, Simplex       == 2, Simplex       == 3)
        print ' TimeSignature : %d %d %d'%(TimeSignature== 1, TimeSignature == 2, TimeSignature == 3)                
        print ' TimeSimplex   : %d %d %d'%(TimeSimplex  == 1, TimeSimplex   == 2, TimeSimplex   == 3)
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        print         
        for i in range(3):
                if ReduceAxes[i] == 1 :
                        print ' %s-axis : positive half stored'%directions[i]
                else:
                        print ' %s-axis : fully         stored'%directions[i]
        print              
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

