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
from src_heph import latex
import sys

# Horizontal line for printing
line  = 80*"-"
#-------------------------------------------------------------------------------
# Hephaestos can be run as
# 
#      python Hephaestos.py #1 #2
#
# where 
#
#   #1 = the name of a functional file, stored in the functionals/ folder.
#   #2 = is a string determining the symmetries of the calculation
#
#-------------------------------------------------------------------------------

heph_name= \
'   =================================================================\n' +\
"   | | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___  | \n" +\
"   | | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __| | \n" +\
"   | |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \ | \n" +\
"   | |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/ | \n" +\
"   |                                                               | \n" +\
"   |  Copyright W.Ryssens & M. Bender                              | \n" +\
"   ================================================================="

print (heph_name)
#-------------------------------------------------------------------------------
# Dealing with the input
inp = sys.argv
if(len(inp) == 1):
    print(" Error: please specify a functional file.")
    print('==============================================================')
    exit()
elif(len(inp) == 2):
    # We select, by default, an EV8-like option
    SYMSTRING = "T,P,Rz,STy"
    REDUCE    = "111"
else:
    SYMSTRING = inp[2]
    REDUCE    = inp[3]

# File containing the definition of the functional
FUNC_FILE  = 'functionals/' + sys.argv[1]

print (line)
print (' Primary input:')
print ('    Functional file: %s'%FUNC_FILE)
print ('    Symmetry string: %s'%SYMSTRING)
print ('    Axis reduction : %s'%REDUCE)
print (line)

#-------------------------------------------------------------------------------
# Path to the original FORTRAN source
SRCPATH = 'src_orig/'
# Path to put the generated source for compilation
GENPATH = 'src/'
#List of FORTRAN files needed for a functional code.
FORTRANFILES=['compilation.f90'  , 'geninfo.f90'      , 'sphericalharmonics.f90',
              'constants.f90'    , 'printing.f90'     , 'HFB.f90','folding.f90' ,
              'nil8.f90'         , 'coulomb.f90'      , 'derivatives.f90'       , 
              'precondition.f90' , 'wavefunctions.f90',
              'hartree-fock.f90' , 'BCS.f90'         , 'pairingcutoffs.f90',
              'temperature_projection.f90', 'momentsofinertia.f90', 
              'particleinabox.f90',
              'densities.f90'    , 'moments.f90'          ,   'pairing.f90'    ,
              'functional.f90'   , 'parameterization.f90' , 'evolution.f90'    , 
              'scfiteration.f90' , 'IO.f90'               , 'tantalus.f90',   
              'transform.f90',  
              'run_single.f90', 'multirun_example.f90', 'timing.f90']

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
# Initialize all of Hephaestos' own modules
heph_symmetries.initsymmetries(SYMSTRING, REDUCE)
print(line)
description = heph_functional.initfunctional(FUNC_FILE)
heph_densities.initdensities()
heph_fields.initfields()

#-------------------------------------------------------------------------------
# Treat all of the source files to a nice dose of preprocessing.
for fname in FORTRANFILES:
     pp.preprocess(fname,SRCPATH,GENPATH)

#-------------------------------------------------------------------------------
# Output all of the relevant things into .tex files.
#print '**************************************************************'
#print ' c)  Hephaestos is creating the .tex files.'
#print '**************************************************************'
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
