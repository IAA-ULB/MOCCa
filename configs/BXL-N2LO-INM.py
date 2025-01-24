# Configuration file for Tantalus compilation
# Functional : N2LO Skyrme-type with BXL-style pairing treatment limited to N2LO terms without external derivatives of densities
FUNC_FILE = 'BXL-N2LO-INM.func'
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
PH_PP_DECOUPL = True
