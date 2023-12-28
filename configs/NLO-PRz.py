# Configuration file for Tantalus compilation
# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
FUNC_FILE = 'NLO.func'
# Symmetries : EV4-style
SYMSTRING = 'T,STy'
REDUCE    = [0,1,0]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,STy'
INREDUCE  = [1,1,0]
# Axis orientation for multipole moments : default
QUANT_AXIS='Z'
SECOND_AXIS=1
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
