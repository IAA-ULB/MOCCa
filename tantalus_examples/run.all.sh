
for script in scripts/run.minimal.sh     scripts/run.pairing.sh    \
              scripts/run.constrained.sh scripts/run.continuing.sh \
              scripts/run.input.sh       scripts/run.blocking.sh   \
              scripts/run.gsk12.sh       scripts/run.gsk12.T.sh    \
              scripts/run.N2LO.sh        scripts/run.cranking.sh   \
              scripts/run.magmoment.sh   scripts/run.Pb208.sh 
do
  echo "Running $script"
  bash $script
done
