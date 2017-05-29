#============================================================================================
#          _______  _______           _______  _______  _______ _________ _______  _______ 
#|\     /|(  ____ \(  ____ )|\     /|(  ___  )(  ____ \(  ____ \\__   __/(  ___  )(  ____ \
#| )   ( || (    \/| (    )|| )   ( || (   ) || (    \/| (    \/   ) (   | (   ) || (    \/
#| (___) || (__    | (____)|| (___) || (___) || (__    | (_____    | |   | |   | || (_____ 
#|  ___  ||  __)   |  _____)|  ___  ||  ___  ||  __)   (_____  )   | |   | |   | |(_____  )
#| (   ) || (      | (      | (   ) || (   ) || (            ) |   | |   | |   | |      ) |
#| )   ( || (____/\| )      | )   ( || )   ( || (____/\/\____) |   | |   | (___) |/\____) |
#|/     \|(_______/|/       |/     \||/     \|(_______/\_______)   )_(   (_______)\_______)
#============================================================================================
# Master script 'forging the chains' of Tantalus. 
#
#
#============================================================================================
from os.path import isfile as isfile
import os
from src_heph import heph_densities, heph_symmetries
from src_heph import preprocess as pp
#-------------------------------------------------------------------------------
# Path to the original FORTRAN source
SRCPATH = 'src_orig/'
# Path to put the generated source for compilation
GENPATH = 'src/'
#List of FORTRAN files to process
FORTRANFILES=['compilation.f90', 'geninfo.f90' , 'constants.f90', 'diag.f90' ,
              'nil8.f90', 'functional.f90',
              'derivatives.f90', 'wavefunctions.f90', 'hartree-fock.f90',
              'densities.f90', 'IO.f90', 'tantalus.f90' ]
#-------------------------------------------------------------------------------
heph_name='============================================================================================\n' + \
'          _______  _______           _______  _______  _______ _________ _______  _______   \n' + \
'|\     /|(  ____ \(  ____ )|\     /|(  ___  )(  ____ \(  ____ \\__   __/(  ___  )(  ____ \  \n' + \
'| )   ( || (    \/| (    )|| )   ( || (   ) || (    \/| (    \/   ) (   | (   ) || (    \/  \n' + \
'| (___) || (__    | (____)|| (___) || (___) || (__    | (_____    | |   | |   | || (_____   \n' + \
'|  ___  ||  __)   |  _____)|  ___  ||  ___  ||  __)   (_____  )   | |   | |   | |(_____  )  \n' + \
'| (   ) || (      | (      | (   ) || (   ) || (            ) |   | |   | |   | |      ) |  \n' + \
'| )   ( || (____/\| )      | )   ( || )   ( || (____/\/\____) |   | |   | (___) |/\____) |  \n' + \
'|/     \|(_______/|/       |/     \||/     \|(_______/\_______)   )_(   (_______)\_______)  \n \n' + \
' Copyright W.Ryssens & M. Bender \n ============================================================================================\n'

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
print ' a)  Hephaestos is initializing its own modules.'
heph_symmetries.initsymmetries()
heph_densities.initdensities()

#-------------------------------------------------------------------------------
# Treat all of the source files to a nice dose of preprocessing.
print ' b)  Hephaestos is processing the template source code.'
for fname in FORTRANFILES:
     print ' * Processing ' + fname
     pp.preprocess(fname,SRCPATH,GENPATH)
    
#-------------------------------------------------------------------------------
# Check if all files got generated correctly. 
print ' c)  Checking that all source code is properly generated.'

FOUND=True    
for fname in FORTRANFILES:
    if(not isfile(GENPATH + fname)):
        print 'You are missing %s%s'%(GENPATH, fname)
        FOUND=False
print "--------------------------------------------------------"        
if(not FOUND):
    print 'Hephaestos did not treat all the source code files.'
else :
    print 'Correct exit. Ready for compilation.'
print "--------------------------------------------------------"


