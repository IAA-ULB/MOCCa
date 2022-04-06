#-------------------------------------------------------------------------------
# Small bash script that generates the public version of Tantalus for the 
# Brussels group with the Brussels EDFs.
#-------------------------------------------------------------------------------

CONFIG=BXL
CONFIG_T=BXL-T
CONFIG_P=BXL-P
CONFIG_TP=BXL-TP
CONFIG_crankX=BXL-crank
CONFIG_T_crankX=BXL-T-crank

SRCPUBLIC=$HOME/Documents/Codes/tantalus_public/src
SRCPUBLIC_T=$HOME/Documents/Codes/tantalus_public/src_T
SRCPUBLIC_P=$HOME/Documents/Codes/tantalus_public/src_P
SRCPUBLIC_TP=$HOME/Documents/Codes/tantalus_public/src_TP

SRCPUBLIC_crankX=$HOME/Documents/Codes/tantalus_public/src_crankX
SRCPUBLIC_T_crankX=$HOME/Documents/Codes/tantalus_public/src_T_crankX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal conserved version; traditional orientation of Qlm

python3 Hephaestos.py $CONFIG &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal conserved version; X-orientation of Qlm

python3 Hephaestos.py $CONFIG_crankX &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_crankX/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal broken version; traditional orientation of Qlm

python3 Hephaestos.py $CONFIG_T &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_T/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Time-reversal broken version; traditional orientation of Qlm

python3 Hephaestos.py $CONFIG_T_crankX &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_T_crankX/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# Parity broken version; traditional orientation of Qlm

python3 Hephaestos.py $CONFIG_P &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_P/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# PT-broken version; traditional orientation of Qlm

python3 Hephaestos.py $CONFIG_TP &>  /dev/null

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRCPUBLIC_TP/

