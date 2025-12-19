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

#-------------------------------------------------------------------------------
# DECLARATIONS 
#-------------------------------------------------------------------------------

# Horizontal line for printing
line  = 80*"-"

heph_name= \
'   ========================================================================== \n' +\
r"   | | | | |   ___   _ __   | |__     __ _    ___   ___  | |_   ___    ___  |"+"\n" +\
r"   | | |_| |  / _ \ | '_ \  | '_ \   / _` |  / _ \ / __| | __| / _ \  / __| |"+"\n" +\
r"   | |  _  | |  __/ | |_) | | | | | | (_| | |  __/ \__ \ | |_ | (_) | \__ \ |"+"\n" +\
r"   | |_| |_|  \___| | .__/  |_| |_|  \__,_|  \___| |___/ \__|  \___/  |___/ |"+"\n" +\
r"   |                |_|                                                     |"+"\n" +\
r"   |                                                                        |"+"\n" +\
r"   |  Copyright W.Ryssens & M. Bender                                       |"+"\n" +\
r"   ==========================================================================="

#List of FORTRAN files needed for a functional code.
FORTRANFILES=['compilation.f90'   , 'geninfo.f90'      , 'sphericalharmonics.f90',
              'constants.f90'     , 'printing.f90'     , 'HFB.f90',
              'HFB_gradient.f90'  , 'HFB_direct.f90'   , 'folding.f90' ,
              'nil8.f90'          , 'coulomb.f90'      , 'derivatives.f90'   ,
              'precondition.f90'  , 'wavefunctions.f90', 'basis_transform.f90',
              'hartree-fock.f90'  , 'BCS.f90'          , 'pairingcutoffs.f90',
              'momentsofinertia.f90', 'hdf5_auxiliary.f90',
              'fission_MOI.f90'   ,  'densities.f90'      ,
              'moments.f90'       , 'pairing.f90'       , 'pairing_strengths.f90',
              'functional.f90'    , 'parameterization.f90' , 'evolution.f90'   ,
              'scfiteration.f90'  , 'IO.f90'               , 'tantalus.f90'    ,
              'transform.f90'     ,  'cranking.f90'        , 'convergence.f90' ,
              'run_single.f90'    , 'multirun_example.f90' , 'timing.f90', 
              'vectors.f90', 'IO_aux.f90', 'IO_wf.f90', 'version.f90', 
              'fam_gmres.f90', 'fam.f90', 'fam_run.f90', 'fam_testing.f90']

#-------------------------------------------------------------------------------
# Actually start Hephaestos
#-------------------------------------------------------------------------------
if(len(sys.argv) != 6):
  
  print ("Running Hephaestos requires specifying either five or four arguments.")
  print (' 1. specifying a configuration file')
  print (' 2. indicate whether compiling a mean-field or FAM executable')
  print (' 3. specifying the chosen option regarding density summation')
  print (' 4. specifying the location of the generated source file')
  print (' 5. name of the source code file to be processed OR "dry-run" to simulate processing of all files')
  print (' Examples:')
  print ('    python Hephaestos.py NLO mf 1 build/NLO/src_mf version.f90 -> will process ONLY version.f90' )
  print ('    python Hephaestos.py NLO mf 1 build/NLO/src_mf dry-run     -> will process ALL .f90 files' )

  sys.exit(1)

config                 = sys.argv[1]
EXETYPE                = sys.argv[2]
DENSITY_SPWF_SUMMATION = int(sys.argv[3])
source_location        = sys.argv[4]
filename               = sys.argv[5]


if(filename != 'dry-run' and not filename.endswith('.f90')):
  print ("The specified filename must end with .f90 or be 'dry-run'.")
  sys.exit(1)

assert(EXETYPE == 'mf' or EXETYPE == 'fam')
if(EXETYPE == 'mf'):  
  fam_active = False
else:
  fam_active = True

if( not os.path.isfile('configs/' + config + '.py')):
  print ("Config file '%s' does not exist."%config)
  sys.exit(1)

if( not os.path.isdir(source_location)):
  print ("The directory you specified for the generated code does not exist.")
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

if(filename == 'dry-run'):
  print (heph_name)
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
SRCPATH = '' #'src_orig/'
GENPATH = source_location 
if(not GENPATH.endswith('/')):
  GENPATH += '/'

#-------------------------------------------------------------------------------
# First we identify all the relevant symmetry options
so    = heph_symmetries.initsymmetries(SYMSTRING,REDUCE,QUANT_AXIS,SECOND_AXIS)
oldso = heph_symmetries.initsymmetries(INSYM,INREDUCE,QUANT_AXIS,SECOND_AXIS)

if(filename == 'dry-run'):
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
description = heph_functional.initfunctional(FUNC_FILE, so, DENSITY_SPWF_SUMMATION, fam_active)
## ... and initialize the fields module
#heph_fields.initfields(so)
#-------------------------------------------------------------------------------
# On to the real business: generating Fortran code.
#for fname in FORTRANFILES:

if( filename == 'dry-run'):
  # This is a dry-run: we will simulate tackling ALL files and print output 
  #    but we will NOT generate any final source code files.
  for fname in FORTRANFILES:
    pp.preprocess(fname,SRCPATH,GENPATH, so, oldso, PH_PP_DECOUPL, fam_active, DENSITY_SPWF_SUMMATION, dry_run=True)
else: 
  # Process the one specific file we have been asked for
  split = filename.split('/')
  fname = split[-1]
  SRCPATH = '/'.join(split[0:-1]) + '/'
  pp.preprocess(fname,SRCPATH,GENPATH, so, oldso, PH_PP_DECOUPL, fam_active, DENSITY_SPWF_SUMMATION, dry_run=False)


sys.exit(0)
