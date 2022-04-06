#-------------------------------------------------------------------------------
# Small bash script that generates the public version of Tantalus for the 
# Lyon group with NLO and N2LO EDFs.
#-------------------------------------------------------------------------------
CONFIG=default
CONFIG_T=NLO-T
CONFIG_P=NLO-P
CONFIG_TP=NLO-TP
CONFIG_crankX=NLO-crank
CONFIG_T_crankX=NLO-T-crank

SRCPUBLIC=$HOME/Documents/Codes/tantalus_lyon/src
SRCPUBLIC_T=$HOME/Documents/Codes/tantalus_lyon/src_T
SRCPUBLIC_P=$HOME/Documents/Codes/tantalus_lyon/src_P
SRCPUBLIC_TP=$HOME/Documents/Codes/tantalus_lyon/src_TP

SRCPUBLIC_crankX=$HOME/Documents/Codes/tantalus_lyon/src_crankX
SRCPUBLIC_T_crankX=$HOME/Documents/Codes/tantalus_lyon/src_T_crankX

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

