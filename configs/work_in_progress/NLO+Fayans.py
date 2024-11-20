# Configuration file for Tantalus compilation
# Functional : NLO Skyrme p-h channel + Fayans pairing
FUNC_FILE = 'NLO+Fayans.func'
# Symmetries : EV8-style
SYMSTRING = 'Rz,T,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='Z'
SECOND_AXIS=1
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = False
