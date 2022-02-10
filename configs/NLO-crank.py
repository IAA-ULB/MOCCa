# Configuration file for Tantalus compilation
# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
# This configuration is different from the default one: we choose a different
# quantisation axis for the multipole moments
FUNC_FILE = 'NLO.func'
# Symmetries : EV8-style
SYMSTRING = 'Rz,T,P,STy'
REDUCE    = [1,1,1]
# Read symmetries: EV8-style
INSYM     = 'Rz,T,P,STy'
INREDUCE  = [1,1,1]
# Axis orientation for multipole moments : default
QUANT_AXIS='X'
SECOND_AXIS=2
# Decouple the particle-hole and particle-particle fields? 
PH_PP_DECOUPL = True
