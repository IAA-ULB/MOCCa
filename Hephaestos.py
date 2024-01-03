#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Master script 'forging the chains' of Tantalus. 
#-------------------------------------------------------------------------------
from os.path  import isfile as isfile
import os
from src_heph import heph_symmetries,heph_densities,heph_functional,heph_fields
from src_heph import preprocess as pp
#from src_heph import latex
import sys, importlib

# Horizontal line for printing
line  = 80*"-"

heph_name= \
'   =================================================================\n' +\
"   | | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___  | \n" +\
"   | | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __| | \n" +\
"   | |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \ | \n" +\
"   | |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/ | \n" +\
"   |              |_|                                              | \n" +\
"   |                                                               | \n" +\
"   |  Copyright W.Ryssens & M. Bender                              | \n" +\
"   ================================================================="

print (heph_name)
#-------------------------------------------------------------------------------
# Dealing with the input:
# a) checking for existence
# b) importing a configuration file
# c) specifying density summation option
#-------------------------------------------------------------------------------
if(len(sys.argv) != 3):
  print ("Running Hephaestos requires specifying two arguments.")
  print (' a) specifying a configuration file')
  print (' b) specifying the chosen option regarding density summation')
  print (' Example:')
  print ('    python Hephaestos.py NLO 1')
  
  sys.exit(1)

config = sys.argv[1]
DENSITY_SPWF_SUMMATION = int(sys.argv[2])

if( not os.path.isfile('configs/' + config + '.py')):
  print ("Config file '%s' does not exist."%config)
  sys.exit(1)

configmod = importlib.import_module('configs.' + config)
# Reassigning the necessary stuff from the configuration file
try:
  FUNC_FILE = 'functionals/' + configmod.FUNC_FILE
except AttributeError:
  print ("Config file does not have a FUNC_FILE attribute.")
  sys.exit(1)

try:
  SYMSTRING = configmod.SYMSTRING
except AttributeError:
  print ("Config file does not have a SYMSTRING attribute.")
  sys.exit(1)

try:
  REDUCE = configmod.REDUCE
except AttributeError:
  print ("Config file does not have a REDUCE attribute.")
  sys.exit(1)

try:
  INSYM = configmod.INSYM
except AttributeError:
  print ("Config file does not have a INSYM attribute.")
  sys.exit(1)

try:
  INREDUCE = configmod.INREDUCE
except AttributeError:
  print ("Config file does not have a INREDUCE attribute.")
  sys.exit(1)

try:
  QUANT_AXIS = configmod.QUANT_AXIS
except AttributeError:
  QUANT_AXIS = 'Z'
  
try:
  SECOND_AXIS = configmod.SECOND_AXIS
except AttributeError:
  SECOND_AXIS = 1

try:
  PH_PP_DECOUPL = configmod.PH_PP_DECOUPL
except AttributeError:
  PH_PP_DECOUPL = True


print (line)
print (' Configuration file                    : %s'%config)
print ('    Functional file                    : %s'%FUNC_FILE)
print ('    PH-PP channel decoupling           : %s'%PH_PP_DECOUPL)
print ('    Symmetry string                    : %s'%SYMSTRING)
print ('    Axis reduction                     : %s'%REDUCE)
print ('    Symmetry string                    : %s'%INSYM)
print ('    Axis reduction                     : %s'%INREDUCE)
print ('    Quantisation axis                  : %s'%QUANT_AXIS)
print ('    Secondary    axis                  : %s'%SECOND_AXIS)
print ('    Sum density derivatives from spwfs : %s'%DENSITY_SPWF_SUMMATION)
print (line)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Path to the original FORTRAN source
SRCPATH = 'src_orig/'
# Path to put the generated source for compilation
GENPATH = 'src/'
#List of FORTRAN files needed for a functional code.
FORTRANFILES=['compilation.f90'   , 'geninfo.f90'      , 'sphericalharmonics.f90',
              'constants.f90'     , 'printing.f90'     , 'HFB.f90',               
              'HFB_gradient.f90'  , 'HFB_direct.f90'   , 'folding.f90' ,
              'nil8.f90'          , 'coulomb.f90'      , 'derivatives.f90'   , 
              'precondition.f90'  , 'wavefunctions.f90', 'basis_transform.f90',
              'hartree-fock.f90'  , 'BCS.f90'          , 'pairingcutoffs.f90',
              'temperature_projection.f90',              'momentsofinertia.f90', 
              'fission_MOI.f90'   , 'particleinabox.f90', 'densities.f90'      , 
              'moments.f90'       , 'pairing.f90'       , 'pairing_strengths.f90',
              'functional.f90'    , 'parameterization.f90' , 'evolution.f90'   , 
              'scfiteration.f90'  , 'IO.f90'               , 'tantalus.f90'    ,   
              'transform.f90'     ,  'cranking.f90'        , 'convergence.f90' , 
              'run_single.f90'    , 'multirun_example.f90' , 'timing.f90']

#-------------------------------------------------------------------------------
# Check for the existence of all the source code files in SRCPATH
FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(SRCPATH + fname)):
        print ('You are missing %s%s'%(SRCPATH, fname))
        FOUND=False
        
if(not FOUND):
    print ('Go find the source files, then come back.')
    quit()

#-------------------------------------------------------------------------------
# First we identify all the relevant symmetry options
so    = heph_symmetries.initsymmetries(SYMSTRING,REDUCE,QUANT_AXIS,SECOND_AXIS)
oldso = heph_symmetries.initsymmetries(INSYM,INREDUCE,QUANT_AXIS,SECOND_AXIS)

print ("  Symmetry information" )
heph_symmetries.printsymmetryoption(so)
print(line)
print ("  Symmetry information" )
heph_symmetries.printsymmetryoption(oldso)
print(line)

# Initialize the densities module, setting up the properties of all the 
# operators
heph_densities.initdensities()
#-------------------------------------------------------------------------------
# Next, we read all the functional information
description = heph_functional.initfunctional(FUNC_FILE, so, DENSITY_SPWF_SUMMATION)
# ... and initialize the fields module
heph_fields.initfields(so)
#-------------------------------------------------------------------------------
# On to the real business: generating Fortran code.
for fname in FORTRANFILES:
#     print ("Preprocessing " + fname)
     pp.preprocess(fname,SRCPATH,GENPATH, so, oldso, PH_PP_DECOUPL, 
                   DENSITY_SPWF_SUMMATION)

#-------------------------------------------------------------------------------
# Output all of the relevant things into .tex files.
#latex.Build(FUNC_FILE, description)

#-------------------------------------------------------------------------------
# Check if all files got generated correctly. 
print (line)
FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(GENPATH + fname)):
        print ('You are missing %s%s'%(GENPATH, fname))
        FOUND=False
if(not FOUND):
    print ('     Hephaestos did not treat all the source code files.')
else :
    print ('     All source files processed.')
    print ('     Ready for compilation.')
    print ('     Happy calculations!')
print (line)
