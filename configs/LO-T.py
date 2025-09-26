# Configuration file for Tantalus compilation
# Functional : LO Skyrme-type (t0t3)
FUNC_FILE = 'LO.func'
# Symmetries : CR8-style
SYMSTRING = 'Rz,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='Z'
SECOND_AXIS=1
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
