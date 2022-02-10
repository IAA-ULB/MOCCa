# Configuration file for Tantalus compilation
# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
FUNC_FILE = 'NLO.func'
# Symmetries : EV8-style
SYMSTRING = 'Rz,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
