# Configuration file for Tantalus compilation
# Functional : N2LO Skyrme-type (think SLy4/5/6) with optional Jmn terms
FUNC_FILE = 'N2LO.func'
# Symmetries : CR4-style
SYMSTRING = 'Rz,STy'
REDUCE    = [1,1,0]
# Read symmetries: EV4-style
INSYM     = 'Rz,T,STy'
INREDUCE  = [1,1,0]
# Axis orientation for multipole moments : default
QUANT_AXIS='Z'
SECOND_AXIS=1
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
