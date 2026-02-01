for f in *.f90
do
    python -m numpy.f2py -c $f -m ${f/.f90/}
done
