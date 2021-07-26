# Small script for building executables for ALL configuration options
# Requires one argument: the compiler.
#
#    build_all.sh gfortran


if [ $# -eq 0 ]; then
    echo "No compiler specified."
fi


for config in configs/*.py
do

c=${config/'configs/'/}
c=${c/'.py'/}

rm -f exec/Tantalus.$c.exe exec/Tantalus.$c.mpi.exe

echo "Compiling configuration $c"
make CONFIG=$c CXX=$1 &> compilation_logs/$c.log
if [ -f "exec/Tantalus.$c.exe" ]; then
 echo "Compilation succesfull."
else
 echo "Compilation failed. Logfile = $c.log"
fi
echo "--------------------------------------------"
done

