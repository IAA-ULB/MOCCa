# Configuration file for Tantalus compilation
# This configuration is different from the BXL one: we choose a different
# quantisation axis for the multipole moments
FUNC_FILE = 'BXL.func'
# Symmetries : EV8-style
SYMSTRING = 'Rz,T,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='X'
SECOND_AXIS=2
