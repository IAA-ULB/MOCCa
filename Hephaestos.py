#===============================================================================
# Master script 'forging the chains' of Tantalus. 
#
#
#===============================================================================
from os.path import isfile as isfile
import os
from src_heph import *

#-------------------------------------------------------------------------------
# Path to the original FORTRAN source
SRCPATH = 'src_orig/'
# Path to put the generated source for compilation
GENPATH = 'src/'
#List of FORTRAN files to process
FORTRANFILES=['compilation.f90', 'geninfo.f90', 'derivatives.f90',
              'wavefunctions.f90', 'tantalus.f90' 
              ]
#-------------------------------------------------------------------------------

print "--------------------------------------------------------"
print 'Welcome to Hephaestos.'
print "Let me generate some code, so that you don't have to."
print "--------------------------------------------------------"

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
# Treat all of the source files to a nice dose of preprocessing.
for fname in FORTRANFILES:
     pp.preprocess(fname,SRCPATH,GENPATH)
    
#-------------------------------------------------------------------------------
# Check if all files got generated correctly. 
FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(GENPATH + fname)):
        print 'You are missing %s%s'%(GENPATH, fname)
        FOUND=False
        
print "--------------------------------------------------------"        
if(not FOUND):
    print 'Hephaestos did not treat all the source code files.'
else :
    print 'Hephaestos did treat all the source code files.'
print "--------------------------------------------------------"


