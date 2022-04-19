# Configuration file for Tantalus compilation
# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
FUNC_FILE = 'NLO.func'
# Symmetries : CR4-style
SYMSTRING = 'Rz,STy'
REDUCE    = [1,1,0]
# Read symmetries: EV4-style
INSYM     = 'Rz,T,STy'
INREDUCE  = [1,1,0]
# Axis orientation for multipole moments : default
QUANT_AXIS='X'
SECOND_AXIS=2
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
