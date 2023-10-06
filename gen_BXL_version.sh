#-------------------------------------------------------------------------------
# Small bash script that generates the public version of Tantalus for the 
# BXL group with NLO and N2LO EDFs.
#-------------------------------------------------------------------------------
for subscript in '' '_T' '_P' '_TP'
do
for EDF       in 'NLO' 'N2LO'
do

EDFSTR="_$EDF"
#echo $EDFSTR$subscript

SRC=$HOME/Documents/Codes/tantalus_public/src$EDFSTR$subscript
echo SRC=$SRC
mkdir -p $SRC

config="BXL-$EDF$subscript$crankX"
config=${config//_/-}
config=${config/'-NLO'/}

#echo $config
ls configs/$config.py

python3 Hephaestos.py $config &> /dev/null
cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRC/

done
done
