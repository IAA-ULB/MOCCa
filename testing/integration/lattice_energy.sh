#------------------------------------------------------------------------------------
# Checking that the difference between Coulomb energy in vacuum and in the lattice
# is equal to the analytical estimate for the spherical nucleus:
# Wl=e^2*Clatt/2 * Z^2/(np*dx) * (1+pi/6/Clatt * rrms^2/(np*dx)^2)
# Clatt=-1.418649 for sc lattice; dx=1.0 for simplicity.
# Ca40 is used for numerical example. Finite nucleon sizes are taken into account.
#
# Order
# 1. The Coulomb energy in vaccum is calculated and compared to the reference value.
#    The charge rms is also checked. BSkG4simpl.param file is used, as it does not 
#    contain any collective corrections.
#
# 2. The Coulomb energy in simple cubic is calculated with the given number of points (10-20 is okey).
#    The correct number of protons is also verified
#
# 3. The difference in two Coulomb energy should be close to the analytical estimate.
#    The error could be due to the perturbation of the Coulomb field induced by the boundary conditions,
#    and due to the mesh discretization. It was checked beforehand that error decreases with the box size.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  Coulomb energy in vacuum              78.527                        50 keV
#  charge rms in vacuum                   3.481                        0.1 fm
#  Lattice energy              <difference of two calcuclations>      100 keV
#  Number of protons                         20                        0.0001
#
#  TODO:
#  - add other observables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#
#  bash lattice_energy.sh -p param -e exec -f exec -n np
#
# -p param: specify a parameterization name (should be without collective corrections)
# -e exec : specify the  suffix of the executable with antiperiodic boundary conditions
#            Example: "BXL.NUCLEI" for "Tantalus.BXL.NUCLEI.exe".
# -f exec : specify the  suffix of the executable with periodic boundary conditions
#            Example: "BXL.PASTA" for "Tantalus.BXL.PASTA.exe".
# -n np   : specify the number of points in the lattice (dx=1.0), so it determines the lattice spacing
#
#  Attention: the exe being called should be able to auto-initialise, i.e. to
#             start from scratch without reading a .wf file!
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : nikolai.shchechilin@ulb.be
# Reference commit hash: e75e6f25c2e2e7223d59043ddf92a6756bb93793
#--------------------------------------------------------------------------------
refEcv=78.527 # Coulomb energy of Ca40 in MeV
refrms=3.4810 # charge rms of Ca40 in fm
refZ=20.0     # number of protons in Ca40 


usage() { echo "Usage: $0 -p param -e exec" 1>&2; exit 1; }

while getopts "p:e:f:n:" opt; do
    case "${opt}" in
        p)
            param=${OPTARG}
            ;;
        e)
            exec_anti=${OPTARG}
            ;;
		f)
            exec_period=${OPTARG}
			;;
		n)
            np=${OPTARG}
			;;
		*)
            usage
            ;;
    esac
done

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A small function to write equivalent input data
write_data()
{
cat << EOF > tant.data
&nucleus
neutrons=20, protons=20
/
&mesh
nx=$2, ny=$2, nz=$2, dx=1.0
/
&func
name_param='$1'
/
&pairing
Type="BCS"
/
&evolution
maxiter=300
/
&scfiteration
/
&wfs
ini_strategy='NILSSON'
nwn = 30, nwp = 30
osc_freq=0.2,0.2,0.2 
/
&IO
InputFilename='init'
Outputfilename='trash'
/
&MomentParam
/
&MomentConstraint
/
&Cranking
/
EOF
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

#----------------------------------------------------------------------------------
# Setting up the reference calculation
setup_test_env "lattice_energy" "$exec_anti" "$param"
write_data  $param 16 
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check1=$?
# get the coulomb energy
Ecv=$(get_coulomb_energy_stdout $outfile)
echo "fast check coul vac" $Ecv $refEcv
# get the rms
rms=$(get_rms_stdout $outfile)
echo "fast check rms" $rms $refrms
# ... and tear down this testing environment.
teardown_test_env
#----------------------------------------------------------------------------------

# Setting up the reference calculation
setup_test_env "lattice_energy" "$exec_period" "$param"
write_data  $param $np
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check2=$?
# get the coulomb energy
Ecl=$(get_coulomb_energy_stdout $outfile)
echo "fast check coul latt" $Ecl
# get Z
protons=$(get_Z_stdout $outfile)
echo "fast check protons" $protons
# ... and tear down this testing environment.
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '------------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------------'
printf ' vacuum calculations   ->  %1d \n' $tantalus_check1
printf ' lattice calculations   ->  %1d \n' $tantalus_check2
# ... and compare with a tolerance of 50 keV to the expected answer
compare_floats $Ecv $refEcv 0.050
check_energy_coul_vac=$?
# ... and compare with a tolerance of 0.1 fm to the expected answer
compare_floats $rms $refrms 0.1
check_rms=$?
# Lattice energy
refEcl=$(echo "1.021402*$protons^2/$np*(1-0.369083*$rms^2/$np^2)" | bc -l )
echo "fast check lattice analyt" $refEcl 
echo "fast check lattice calc:"
echo "$Ecv-$Ecl" |bc -l 
# ... and compare the difference with a tolerance of 100 keV to lattice energy
compare_floats $Ecv-$Ecl $refEcl 0.1
check_energy_coul_lat=$?
# ... and compare with a tolerance of 0.001 to the expected answer
compare_floats $protons $refZ 0.0001
check_Z=$?

echo '------------------------------------------------'
printf ' Coulomb energy in vacuum consistency->  %1d \n' $check_energy_coul_vac
echo '------------------------------------------------'
echo '------------------------------------------------'
printf ' Charge rms consistency              ->  %1d \n' $check_rms
echo '------------------------------------------------'
echo '------------------------------------------------'
printf ' Lattice energy consistency          ->  %1d \n' $check_energy_coul_lat
echo '------------------------------------------------'
echo '------------------------------------------------'
printf ' Number of protons consistency       ->  %1d \n' $check_Z
echo '------------------------------------------------'

exitcode=$(( $check_energy_coul_vac || $check_energy_coul_lat || $check_rms || $check_Z ))
echo '------------------------------------------------'
printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
