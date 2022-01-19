
for script in scripts/run.minimal.sh     scripts/run.pairing.sh    \
              scripts/run.constrained.sh scripts/run.continuing.sh \
              scripts/run.input.sh       scripts/run.blocking.sh   \
              scripts/run.gradient.sh    scripts/run.cranking.sh   \
              scripts/run.BSkG.sh        scripts/run.BSkG.T.sh     \
              scripts/run.N2LO.sh        scripts/run.cranking.sh   \
              scripts/run.magmoment.sh   scripts/run.P.sh          \
              scripts/run.Pb208.sh 
do
  echo "Running $script"
  bash $script
done
