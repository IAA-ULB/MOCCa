# Small script for building executables for ALL configuration options
# Requires one argument: the compiler.
#
#    build_all.sh gfortran

# Preparation
if [ $# -eq 0 ]; then
    echo "No compiler specified."
    exit
fi

mkdir -p compilation_logs

# Loop over all configurations
for config in configs/*.py
do

c=${config/'configs/'/}
c=${c/'.py'/}

# remove older executables
rm -f exec/Tantalus.$c.exe exec/Tantalus.$c.mpi.exe

echo "Compiling configuration $c"
make CONFIG=$c CXX=$1 &> compilation_logs/$c.log
if [ -f "exec/Tantalus.$c.exe" ]; then
 echo "Compilation succesful."
else
 echo "Compilation failed. Logfile = $c.log"
fi
echo "--------------------------------------------"
done

