# Small bash script that generates the public version of Tantalus
CONFIG=DD-switch
CONFIG_T=DD-switch-T
CONFIG_P=DD-switch-P
CONFIG_TP=DD-switch-TP

SRCPUBLIC=$HOME/Documents/Codes/tantalus_public/src
SRCPUBLIC_T=$HOME/Documents/Codes/tantalus_public/src_T
SRCPUBLIC_P=$HOME/Documents/Codes/tantalus_public/src_P
SRCPUBLIC_TP=$HOME/Documents/Codes/tantalus_public/src_TP

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal conserved version

python3 Hephaestos.py $CONFIG

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal broken version

python3 Hephaestos.py $CONFIG_T

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_T/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Parity broken version

python3 Hephaestos.py $CONFIG_P

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_P/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# PT-broken version

python3 Hephaestos.py $CONFIG_TP

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_TP/

