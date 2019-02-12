#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module that determines the symmetries imposed by Hephaestos on 
# Tantalus/Prometheus.
#
#
#
#
#
#-------------------------------------------------------------------------------

#===============================================================================
#
# Whether or not signature is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
Signature     = 0
#
# Whether or not parity is a quantum number of the spwfs. 
#   0 Broken
#   1 Conserved
Parity        = 0
#
# Whether or not time-simplex is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
#
TimeSimplex   = 0
#
# Whether or not time-signature is a quantum number of the spwfs, 
# and which Cartesian direction it is applied on.
#   0 Broken
#   1 Conserved, X
#   2 Conserved, Y
#   3 Conserved, Z
#
TimeSignature = 0
#
# Whether or not time-reversal is conserved. Can be 0 (broken) or 1 (conserved).
TimeReversal  = 1
#
#===============================================================================
# Whether you want a reduction of each of the Cartesian axes. 
# Note that not all choices are compatible with all symmetry choices. 
# Hephaestos will complain if it feels something is amiss.
# Reduce_x/y/z = 1 if reduced, 0 if fully stored. 
ReduceAxes= [1,1,1]

def initsymmetries():   
        global Signature, Simplex, Parity, TimeSimplex, TimeReversal
        global ReduceAxes     

        # EV8 symmetries for now
        Signature = 3    ; Simplex = 0 ; Parity = 1 ; TimeSimplex = 2;
        TimeReversal = 1 ; TimeSignature = 0

        CheckSymmetryChoice()
        CheckReduction()
        PrintSymmetry()

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

def CheckSymmetryChoice() :
    #---------------------------------------------------------------------------
    # Checking whether the choices made for the symmetries are valid.
    global Signature, Simplex, Parity, TimeSimplex, TimeReversal, TimeSignature
    
    # First off, does anything take illegal values 
    if(Signature > 3 or Signature < 0 ):
            print 'Illegal value for signature.'
            exit()
    if(Simplex > 3 or Simplex < 0 ):
            print 'Illegal value for simplex.'
            exit()
    if(TimeSimplex > 3 or TimeSimplex < 0 ):
            print 'Illegal value for time-simplex.'
            exit()
    if(TimeSignature > 3 or TimeSignature < 0 ):
            print 'Illegal value for time-signature.'
            exit()
    if(Parity > 1 or Parity < 0 ):
            print 'Illegal value for parity.'
            exit()
    if(TimeReversal > 1 or TimeReversal < 0 ):
            print 'Illegal value for time-reversal.'
            exit()

def CheckReduction():
    #---------------------------------------------------------------------------
    # Check whether the asked-for reduction of axes is valid, given the 
    # particular combination of symmetries asked for.

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
     print 'You asked for more reduction of axes than the symmetries can cover.'
     exit()
    elif(s < p):
     print 'The symmetries you demanded provide more reduction of axes ' \
         + 'than you asked for.'
     exit()

    #-----------------------------------------------------------------------
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
