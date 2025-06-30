# Configuration file for Tantalus compilation
# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
FUNC_FILE = 'BXL.func'
# Symmetries : EV4-style
SYMSTRING = 'Rz,T,STy'
REDUCE    = [1,1,0]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='X'
SECOND_AXIS=2
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
