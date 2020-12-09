


for script in scripts/run.minimal.sh     scripts/run.pairing.sh    \
              scripts/run.constrained.sh scripts/run.continuing.sh \
              scripts/run.gsk12.sh       scripts/run.Pb208.sh      \
              scripts/run.N2LO.sh
do
  echo "Running $script"
  bash $script
done
