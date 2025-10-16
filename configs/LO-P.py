# Configuration file for Tantalus compilation
# Functional : LO Skyrme-type
FUNC_FILE = 'LO.func'
# Symmetries : EV4-style
SYMSTRING = 'Rz,T,STy'
REDUCE    = [1,1,0]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='Z'
SECOND_AXIS=1
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
