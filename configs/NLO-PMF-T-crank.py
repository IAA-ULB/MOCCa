# Configuration file for Tantalus compilation
# Multipole moments quantised along the X direction
FUNC_FILE = 'NLO.func'
# Symmetries : EV8-style
SYMSTRING = 'Rz,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : X
QUANT_AXIS='X'
SECOND_AXIS=2
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = False
