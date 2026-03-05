#python -m numpy.f2py -c pow.f90 -m pow

for f in sp*.f90
do
    python -m numpy.f2py -c $f -m ${f/.f90/}
done
