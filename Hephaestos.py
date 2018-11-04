#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Master script 'forging the chains' of Tantalus. 
#
#
#-------------------------------------------------------------------------------
from os.path import isfile as isfile
import os
from src_heph import heph_densities,heph_symmetries,heph_functional,heph_fields
from src_heph import preprocess as pp
from src_heph import latex
import sys


#-------------------------------------------------------------------------------
# Path to the original FORTRAN source
SRCPATH = 'src_orig/'
# Path to put the generated source for compilation
GENPATH = 'src/'
#List of FORTRAN files needed for a functional code.
FORTRANFILES=['compilation.f90'  , 'geninfo.f90' , 'sphericalharmonics.f90',
              'constants.f90'   , 'printing.f90',      'HFB.f90','folding.f90',
              'diag.f90'         , 'nil8.f90'         , 'coulomb.f90'     ,
              'derivatives.f90'  , 'precondition.f90' , 'wavefunctions.f90',
              'hartree-fock.f90' , 'BCS.f90' , 'pairingcutoffs.f90',
              'densities.f90'    , 'moments.f90'          ,   'pairing.f90'    ,
              'functional.f90'   , 'parameterization.f90' , 'evolution.f90'    , 
              'scfiteration.f90' , 'IO.f90'           , 'tantalus.f90' ]

# File containing the definition of the functional
FUNC_FILE      = 'functionals/' + sys.argv[1]

#-------------------------------------------------------------------------------
heph_name= \
'==============================================================\n' +\
"| | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___  \n" +\
"| |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __| \n" +\
"|  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \ \n" +\
"|_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/ \n" +\
" Copyright W.Ryssens & M. Bender \n"                              +\
"==============================================================\n"

print heph_name

#-------------------------------------------------------------------------------
# Check for the existence of all the source code files in SRCPATH
FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(SRCPATH + fname)):
        print 'You are missing %s%s'%(SRCPATH, fname)
        FOUND=False
        
if(not FOUND):
    print 'Go find the source files, then come back.'
    quit()

#-------------------------------------------------------------------------------
# Initialize all of Hephaestos' own modules
print '**************************************************************'
print ' a)  Hephaestos is initializing its own modules.'
heph_symmetries.initsymmetries()
heph_functional.initfunctional(FUNC_FILE)
heph_densities.initdensities()
heph_fields.initfields()
print '**************************************************************'
#-------------------------------------------------------------------------------
# Output all of the relevant things into .tex files.

latex.FillDensities()
exit()
#-------------------------------------------------------------------------------
# Treat all of the source files to a nice dose of preprocessing.
print 
print '**************************************************************'
print ' b)  Hephaestos is processing the template source code.'
print '**************************************************************'
for fname in FORTRANFILES:
     print ' * Processing ' + fname
     pp.preprocess(fname,SRCPATH,GENPATH)
    
#-------------------------------------------------------------------------------
# Check if all files got generated correctly. 
print 
print '**************************************************************'
print ' c)  Checking that all source code is properly generated.'

FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(GENPATH + fname)):
        print 'You are missing %s%s'%(GENPATH, fname)
        FOUND=False
if(not FOUND):
    print '     Hephaestos did not treat all the source code files.'
else :
    print '     Correct exit. Ready for compilation.'
    print '     Happy optimizing!'
print '**************************************************************'
