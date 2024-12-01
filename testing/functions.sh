################################################################################
# Boilerplate functions for writing a Tantalus testing script
#
################################################################################

# Hard-coded variables that guarantee the correct working of the functions below
# TODO: find some way to not hardcode this
EXECDIR=$HOME/Documents/Codes/Tantalus/exec/
PARAMDIR=$HOME/Documents/Codes/Tantalus/parameterizations/

setup_test_env () {
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setting up for running a test: 
# 1. define a few standard environment variables 
# 2. create a work directory 
# 3. copy the relevant executable and .param file there
#
# Arguments are:
#  $1 => configuration file name, or rather the X in Tantalus.X.exe
#  $2 => parameterization name, or rather the X in X.param
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#1. environment variables
exe="Tantalus.$1.exe"  # full name of the executable
param="$2"

#2. create working directory
if [ ! -d "work/" ]; then
  mkdir work
fi

#3. copy executable and parameterization file
cp $EXECDIR/$exe             work/
cp $PARAMDIR/"$param.param"  work/

cd work
}

tantalus_error_codes () {
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Act on the exit code returned by Tantalus
# TODO: do something meaningful with these error codes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $1 -eq 0 ]; then
  sleep 0
  #echo "Tantalus ran succesfully."
fi
}

get_total_energy_stdout (){
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get the final total energy from the Tantalus STDOUT as a float
#
# Input:
#    $1: filename of tantalus STDOUT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  echo `grep "Total energy" $1  | tail -1 | grep -oE '[+-][0-9]+([.][0-9]+)?'`
}

compare_floats (){
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare two floating point numbers for equality within a given tolerance.
# I.e. the exit status of this function is whether or not the following is true
#    | $1 - $2 |  <  $3 
# Caveat: this function relies on the limited bash basic calculator.
#
# Input:
#   $1, $2 : numbers to be compared
#   $3     : tolerance for the comparison
#
# Return codes:
#   0      : comparison is true
#   1      : comparison is false
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Calculate the difference between both floats with the basic calculator
difference=$(echo "$1 - $2" | bc )
# Take the absolute value =  remove the first occurence of "-"
difference=${difference#-}

# Use the bc calculator again to compare the difference to a tolerance
if [ 1 -eq "$(echo "$difference < $3 " | bc)" ]
then 
 return 0 # explicit returns instead of exits to hand control back
else
 return 1
fi
}
