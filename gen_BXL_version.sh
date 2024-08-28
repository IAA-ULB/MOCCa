#-------------------------------------------------------------------------------
# Small bash script that generates the public version of Tantalus for the 
# BXL group with NLO and N2LO EDFs.
#-------------------------------------------------------------------------------
#set -v
for subscript in '' '_T' '_P' '_TP'
do
for EDF       in 'NLO' 'N2LO'
do

EDFSTR="_$EDF"
#echo $EDFSTR$subscript

SRC=$HOME/Documents/Codes/tantalus_public/src$EDFSTR$subscript
SRC=$SRC
mkdir -p $SRC

config="BXL-$EDF$subscript$crankX"
config=${config//_/-}
config=${config/'-NLO'/}

printf "%-30s ->         %-30s \n" configs/$config.py $SRC

# Do the Hephaestos work silently, specifying CONFIG and DENSUM options, ...

python3 Hephaestos.py $config 0 &> /dev/null
# ... but check its return code
if [ $? -eq 1 ]; then
   echo $?
   echo "Hephaestos exited unsuccesfully."
   exit
fi

cp src/tantalus.f90 src/tantalus.version.f90
sed -i.bak 's/VERSION1/"${GIT_INFO1}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION2/"${GIT_INFO2}"/' src/tantalus.version.f90 
sed -i.bak 's/VERSION3/"${GIT_INFO3}"/' src/tantalus.version.f90 
rm src/tantalus.version.f90.bak

cp src/*.f90 $SRC/

done
done
